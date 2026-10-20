import sys, time, importlib.util, torch, numpy as np
W = 'source/weights/model/'
spec = importlib.util.spec_from_file_location('sei', W + 'sei.py'); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
net = mod.Sei(); sd = torch.load(W + 'sei.pth', map_location='cpu')
sd = {k.replace('module.model.', '').replace('module.', ''): v for k, v in sd.items()}
print(net.load_state_dict(sd, strict=False))
net.eval()
with torch.no_grad(): net(torch.zeros(1, 4, 4096))  # builds the lazy spline matrix on cpu
for dev in sys.argv[1:]:
    m = net.to(dev)
    for mod_ in m.modules():
        if getattr(mod_, '_spline_tr', None) is not None: mod_._spline_tr = mod_._spline_tr.to(dev)
    x = torch.eye(4)[torch.randint(0, 4, (32, 4096))].permute(0, 2, 1).contiguous().to(dev)
    with torch.no_grad():
        m(x[:4]); t = time.time(); y = m(x); getattr(torch, dev).synchronize() if dev == 'mps' else None
    print(dev, 'batch 32:', round(time.time() - t, 2), 's', tuple(y.shape), float(y.mean()))
