"""read the first sheet of an xlsx whose zip central directory is damaged (GEO GSE152404), via local file headers; writes a tsv."""
import struct, sys, zlib, xml.etree.ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def members(data):
    i = 0
    while data[i:i + 4] == b'PK\x03\x04':
        _, flag, comp, _, _, _, cs, _, nl, el = struct.unpack('<HHHHHIIIHH', data[i + 4:i + 30])
        name, start = data[i + 30:i + 30 + nl].decode(), i + 30 + nl + el
        assert comp == 8 and not flag & 8, 'unexpected zip layout'
        yield name, zlib.decompress(data[start:start + cs], -15)
        i = start + cs

def col(ref):
    n = 0
    for ch in ref:
        if not ch.isalpha(): break
        n = n * 26 + ord(ch) - 64
    return n - 1

src, dst = sys.argv[1:3]
z = dict(members(open(src, 'rb').read()))
ss = [''.join(t.itertext()) for t in ET.fromstring(z['xl/sharedStrings.xml']).findall('m:si', NS)]
with open(dst, 'w') as fh:
    for r in ET.fromstring(z['xl/worksheets/sheet1.xml']).find('m:sheetData', NS).findall('m:row', NS):
        cells = {}
        for c in r.findall('m:c', NS):
            v = c.find('m:v', NS)
            cells[col(c.get('r'))] = '' if v is None else (ss[int(v.text)] if c.get('t') == 's' else v.text)
        fh.write('\t'.join(cells.get(j, '') for j in range(max(cells) + 1 if cells else 0)) + '\n')
