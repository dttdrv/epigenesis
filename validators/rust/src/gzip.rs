//! Pure-std RFC 1952 (gzip) + RFC 1951 (DEFLATE) decompressor with CRC-32 and size bounds.

const MAX_DECOMPRESSED_BYTES: usize = 64 * 1024 * 1024;

fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xffff_ffff;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xedb8_8320 & mask);
        }
    }
    !crc
}

struct BitReader<'a> {
    data: &'a [u8],
    byte_pos: usize,
    bit_buf: u32,
    bit_cnt: u8,
}

impl<'a> BitReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            byte_pos: 0,
            bit_buf: 0,
            bit_cnt: 0,
        }
    }

    fn ensure(&mut self, n: u8) -> Result<(), String> {
        while self.bit_cnt < n {
            if self.byte_pos >= self.data.len() {
                return Err("unexpected end of deflate stream".to_string());
            }
            self.bit_buf |= (self.data[self.byte_pos] as u32) << self.bit_cnt;
            self.byte_pos += 1;
            self.bit_cnt += 8;
        }
        Ok(())
    }

    fn read_bits(&mut self, n: u8) -> Result<u32, String> {
        if n == 0 {
            return Ok(0);
        }
        self.ensure(n)?;
        let val = self.bit_buf & ((1u32 << n) - 1);
        self.bit_buf >>= n;
        self.bit_cnt -= n;
        Ok(val)
    }

    fn align_byte(&mut self) {
        let drop = self.bit_cnt % 8;
        self.bit_buf >>= drop;
        self.bit_cnt -= drop;
    }

    fn read_aligned_u16(&mut self) -> Result<u16, String> {
        self.align_byte();
        let lo = self.read_bits(8)?;
        let hi = self.read_bits(8)?;
        Ok(((hi << 8) | lo) as u16)
    }
}

struct HuffTable {
    counts: [u16; 16],
    symbols: Vec<u16>,
}

impl HuffTable {
    fn from_lengths(lengths: &[u8]) -> Result<Self, String> {
        let mut counts = [0u16; 16];
        for &len in lengths {
            if len as usize >= 16 {
                return Err("Huffman code length > 15".to_string());
            }
            if len > 0 {
                counts[len as usize] += 1;
            }
        }
        let mut left = 1i32;
        for bits in 1..=15 {
            left <<= 1;
            left -= counts[bits] as i32;
            if left < 0 {
                return Err("over-subscribed Huffman tree".to_string());
            }
        }
        let mut offsets = [0u16; 16];
        for bits in 1..15 {
            offsets[bits + 1] = offsets[bits] + counts[bits];
        }
        let total: usize = counts.iter().map(|&c| c as usize).sum();
        let mut symbols = vec![0u16; total];
        for (sym, &len) in lengths.iter().enumerate() {
            if len > 0 {
                let idx = offsets[len as usize] as usize;
                offsets[len as usize] += 1;
                symbols[idx] = sym as u16;
            }
        }
        Ok(Self { counts, symbols })
    }

    fn decode(&self, br: &mut BitReader<'_>) -> Result<u16, String> {
        let mut code: i32 = 0;
        let mut first: i32 = 0;
        let mut index: i32 = 0;
        for len in 1..=15 {
            code |= br.read_bits(1)? as i32;
            let count = self.counts[len] as i32;
            if code - count < first {
                return Ok(self.symbols[(index + (code - first)) as usize]);
            }
            index += count;
            first += count;
            first <<= 1;
            code <<= 1;
        }
        Err("invalid Huffman code in deflate stream".to_string())
    }
}

fn fixed_tables() -> (HuffTable, HuffTable) {
    let mut lit_len = [0u8; 288];
    for item in lit_len.iter_mut().take(144) {
        *item = 8;
    }
    for item in lit_len.iter_mut().take(256).skip(144) {
        *item = 9;
    }
    for item in lit_len.iter_mut().take(280).skip(256) {
        *item = 7;
    }
    for item in lit_len.iter_mut().take(288).skip(280) {
        *item = 8;
    }
    let dist = [5u8; 32];
    (
        HuffTable::from_lengths(&lit_len).expect("valid fixed litlen"),
        HuffTable::from_lengths(&dist).expect("valid fixed dist"),
    )
}

fn dynamic_tables(br: &mut BitReader<'_>) -> Result<(HuffTable, HuffTable), String> {
    const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
    let hlit = br.read_bits(5)? as usize + 257;
    let hdist = br.read_bits(5)? as usize + 1;
    let hclen = br.read_bits(4)? as usize + 4;
    let mut cl_lens = [0u8; 19];
    for i in 0..hclen {
        cl_lens[ORDER[i]] = br.read_bits(3)? as u8;
    }
    let cl_table = HuffTable::from_lengths(&cl_lens)?;
    let total = hlit + hdist;
    let mut lens = Vec::with_capacity(total);
    while lens.len() < total {
        let sym = cl_table.decode(br)?;
        match sym {
            0..=15 => lens.push(sym as u8),
            16 => {
                let prev = *lens.last().ok_or_else(|| "repeat with no previous length".to_string())?;
                let rep = br.read_bits(2)? as usize + 3;
                for _ in 0..rep {
                    lens.push(prev);
                }
            }
            17 => {
                let rep = br.read_bits(3)? as usize + 3;
                for _ in 0..rep {
                    lens.push(0);
                }
            }
            18 => {
                let rep = br.read_bits(7)? as usize + 11;
                for _ in 0..rep {
                    lens.push(0);
                }
            }
            _ => return Err("invalid code length symbol".to_string()),
        }
    }
    if lens.len() != total {
        return Err("code length overflow in dynamic Huffman header".to_string());
    }
    Ok((
        HuffTable::from_lengths(&lens[..hlit])?,
        HuffTable::from_lengths(&lens[hlit..])?,
    ))
}

fn inflate_stream(stream: &[u8], out: &mut Vec<u8>) -> Result<usize, String> {
    const LEN_BASE: [u16; 29] = [
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115,
        131, 163, 195, 227, 258,
    ];
    const LEN_EXTRA: [u8; 29] = [
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
    ];
    const DIST_BASE: [u16; 30] = [
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
        2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
    ];
    const DIST_EXTRA: [u8; 30] = [
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12,
        13, 13,
    ];

    let mut br = BitReader::new(stream);
    loop {
        let bfinal = br.read_bits(1)?;
        let btype = br.read_bits(2)?;
        if btype == 0 {
            let len = br.read_aligned_u16()?;
            let nlen = br.read_aligned_u16()?;
            if len != !nlen {
                return Err("corrupt stored block length in deflate".to_string());
            }
            if out.len() + len as usize > MAX_DECOMPRESSED_BYTES {
                return Err("decompressed gzip payload exceeds limit".to_string());
            }
            for _ in 0..len {
                out.push(br.read_bits(8)? as u8);
            }
        } else if btype == 1 || btype == 2 {
            let (lit_tab, dist_tab) = if btype == 1 {
                fixed_tables()
            } else {
                dynamic_tables(&mut br)?
            };
            loop {
                let sym = lit_tab.decode(&mut br)?;
                if sym < 256 {
                    if out.len() + 1 > MAX_DECOMPRESSED_BYTES {
                        return Err("decompressed gzip payload exceeds limit".to_string());
                    }
                    out.push(sym as u8);
                } else if sym == 256 {
                    break;
                } else {
                    let li = (sym - 257) as usize;
                    if li >= LEN_BASE.len() {
                        return Err("invalid length symbol in deflate".to_string());
                    }
                    let length = LEN_BASE[li] as usize + br.read_bits(LEN_EXTRA[li])? as usize;
                    let dsym = dist_tab.decode(&mut br)? as usize;
                    if dsym >= DIST_BASE.len() {
                        return Err("invalid distance symbol in deflate".to_string());
                    }
                    let dist = DIST_BASE[dsym] as usize + br.read_bits(DIST_EXTRA[dsym])? as usize;
                    if dist == 0 || dist > out.len() {
                        return Err("invalid back-reference distance in deflate".to_string());
                    }
                    if out.len() + length > MAX_DECOMPRESSED_BYTES {
                        return Err("decompressed gzip payload exceeds limit".to_string());
                    }
                    for _ in 0..length {
                        let b = out[out.len() - dist];
                        out.push(b);
                    }
                }
            }
        } else {
            return Err("reserved deflate block type".to_string());
        }
        if bfinal == 1 {
            break;
        }
    }
    br.align_byte();
    let unused_bytes = (br.bit_cnt / 8) as usize;
    Ok(br.byte_pos - unused_bytes)
}

pub fn decompress_gzip(raw: &[u8]) -> Result<Vec<u8>, String> {
    let mut pos = 0;
    let mut out = Vec::new();
    while pos < raw.len() {
        if pos + 10 > raw.len() || raw[pos] != 0x1f || raw[pos + 1] != 0x8b || raw[pos + 2] != 0x08 {
            return Err("invalid gzip member header".to_string());
        }
        let flg = raw[pos + 3];
        if (flg & 0xe0) != 0 {
            return Err("unsupported reserved gzip header flags".to_string());
        }
        pos += 10;
        if (flg & 0x04) != 0 {
            if pos + 2 > raw.len() {
                return Err("truncated gzip FEXTRA header".to_string());
            }
            let xlen = u16::from_le_bytes([raw[pos], raw[pos + 1]]) as usize;
            pos += 2 + xlen;
        }
        if (flg & 0x08) != 0 {
            while pos < raw.len() && raw[pos] != 0 {
                pos += 1;
            }
            pos += 1;
        }
        if (flg & 0x10) != 0 {
            while pos < raw.len() && raw[pos] != 0 {
                pos += 1;
            }
            pos += 1;
        }
        if (flg & 0x02) != 0 {
            pos += 2;
        }
        if pos > raw.len() {
            return Err("truncated gzip header".to_string());
        }
        let start_out = out.len();
        let consumed = inflate_stream(&raw[pos..], &mut out)?;
        pos += consumed;
        if pos + 8 > raw.len() {
            return Err("truncated gzip trailer".to_string());
        }
        let expected_crc = u32::from_le_bytes([raw[pos], raw[pos + 1], raw[pos + 2], raw[pos + 3]]);
        let expected_isize = u32::from_le_bytes([raw[pos + 4], raw[pos + 5], raw[pos + 6], raw[pos + 7]]);
        pos += 8;
        let member_bytes = &out[start_out..];
        if crc32(member_bytes) != expected_crc {
            return Err("gzip CRC-32 mismatch".to_string());
        }
        if (member_bytes.len() as u32) != expected_isize {
            return Err("gzip ISIZE mismatch".to_string());
        }
    }
    if out.is_empty() && raw.is_empty() {
        return Err("empty gzip input".to_string());
    }
    Ok(out)
}
