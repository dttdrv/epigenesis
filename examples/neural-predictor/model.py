"""One position-aware block with a scalar regulatory-score interface."""

import torch
from torch import nn
from torch.nn import functional as F


def rational_gate(a, b):
    scale = torch.maximum(torch.ones_like(a), torch.maximum(a.abs(), b.abs())).detach()
    p, q, t = a/scale, b/scale, scale.reciprocal()
    return p*(t+q)/(t*t+p*p+q*q)


class RegulatoryAttention(nn.Module):
    def __init__(self, *, length, width, heads, hidden, kernel):
        super().__init__()
        if any(type(v) is not int or v <= 0 for v in (length, width, heads, hidden, kernel)):
            raise ValueError('model dimensions must be positive integers')
        if width % heads or kernel % 2 != 1 or kernel > length:
            raise ValueError('heads must divide width and the odd kernel must fit the sequence')
        self.length, self.width, self.heads = length, width, heads
        position = torch.arange(length)
        self.register_buffer('distance', position[None, :]-position[:, None], persistent=False)
        self.stem = nn.Conv1d(4, width, kernel, padding=kernel//2, bias=False)
        self.attention_norm = nn.LayerNorm(width)
        self.qkv = nn.Linear(width, 3*width, bias=False)
        self.project = nn.Linear(width, width, bias=False)
        self.position_weight = nn.Parameter(torch.zeros(heads, 2))
        self.feed_norm = nn.LayerNorm(width)
        self.feed_up = nn.Linear(width, 2*hidden, bias=False)
        self.feed_down = nn.Linear(hidden, width, bias=False)
        self.final_norm = nn.LayerNorm(width)
        self.readout = nn.Linear(width, 1, bias=False)

    def attend(self, value):
        q, k, v = self.qkv(value).reshape(
            len(value), self.length, 3, self.heads, self.width//self.heads
        ).permute(2, 0, 3, 1, 4).unbind(0)
        distance = self.distance.to(value.dtype)
        bias = (self.position_weight[:, 0, None, None]*distance
                + self.position_weight[:, 1, None, None]*distance.abs())
        attended = F.scaled_dot_product_attention(q, k, v, attn_mask=bias[None], dropout_p=0.0)
        return self.project(attended.transpose(1, 2).reshape(len(value), self.length, self.width))

    def forward(self, value):
        if value.ndim != 3 or not value.shape[0] or value.shape[1:] != (4, self.length):
            raise ValueError('expected a nonempty batch of four-channel full-length sequences')
        value = F.silu(self.stem(value)).transpose(1, 2)
        value = value+self.attend(self.attention_norm(value))
        a, b = self.feed_up(self.feed_norm(value)).chunk(2, dim=-1)
        value = value+self.feed_down(rational_gate(a, b))
        return self.readout(self.final_norm(value)).mean(dim=1).squeeze(-1)
