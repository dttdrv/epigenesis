"""Replay a measured neural perturbation comparison from original source bytes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from brainc._io import MAX_DECOMPRESSED_BYTES, load_json_object, read_regular_file

HERE = Path(__file__).resolve().parent
IMPLEMENTATION_SHA256 = hashlib.sha256(read_regular_file(HERE/'evidence.py')).hexdigest()
CASE_SHA256 = 'b35cbcce124ed3e1a159a6b4c3af1dc51eff47eaf42591e3342716a2b28e3ecb'
X = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
R = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n').encode()


def load_case(base: Path = HERE) -> tuple[dict, dict[str, bytes]]:
    raw = read_regular_file(base/'case.json')
    if sha(raw) != CASE_SHA256:
        raise ValueError('reviewed case identity mismatch')
    case = load_json_object(raw, 'case')
    sources = {}
    for source in case['sources']:
        name = source['file']
        if Path(name).name != name or name in sources:
            raise ValueError('invalid source filename')
        content = read_regular_file(base/'data'/name, label=name)
        if sha(content) != source['sha256']:
            raise ValueError(f'source identity mismatch: {name}')
        sources[name] = content
    return case, sources


def office(raw: bytes) -> dict[str, bytes]:
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError('duplicate Office archive members')
        if sum(i.file_size for i in archive.infolist()) > MAX_DECOMPRESSED_BYTES:
            raise ValueError('Office archive exceeds decompressed byte limit')
        return {name: archive.read(name) for name in names}


def xml(raw: bytes) -> ET.Element:
    if b'<!DOCTYPE' in raw.upper() or b'<!ENTITY' in raw.upper():
        raise ValueError('XML declarations are unsupported')
    return ET.fromstring(raw)


def workbook(raw: bytes) -> dict[str, dict]:
    parts = office(raw)
    strings = []
    if 'xl/sharedStrings.xml' in parts:
        strings = [''.join(t.text or '' for t in s.iter(X+'t'))
                   for s in xml(parts['xl/sharedStrings.xml'])]
    relations = {r.attrib['Id']: r.attrib['Target']
                 for r in xml(parts['xl/_rels/workbook.xml.rels'])}
    result = {}
    for sheet in xml(parts['xl/workbook.xml']).iter(X+'sheet'):
        name = sheet.attrib['name']
        if name in result:
            raise ValueError('duplicate worksheet name')
        target = relations[sheet.attrib[R+'id']]
        path = PurePosixPath(target.lstrip('/') if target.startswith('/') else 'xl/'+target)
        if '..' in path.parts:
            raise ValueError('unsupported worksheet relationship')
        cells = {}
        for cell in xml(parts[str(path)]).iter(X+'c'):
            address = cell.attrib['r']
            if re.fullmatch(r'[A-Z]+[1-9][0-9]*', address) is None or address in cells:
                raise ValueError('invalid or duplicate cell address')
            if cell.find(X+'f') is not None:
                raise ValueError('formula cells require evaluation and are unsupported')
            kind = cell.get('t', 'n')
            value = cell.findtext(X+'v')
            if kind == 'inlineStr':
                parsed = ''.join(t.text or '' for t in cell.iter(X+'t'))
            elif value is None:
                parsed = None
            elif kind == 's':
                index = int(value)
                if not 0 <= index < len(strings):
                    raise ValueError('invalid shared string index')
                parsed = strings[index]
            elif kind == 'n':
                parsed = int(value) if re.fullmatch(r'-?[0-9]+', value) else float(value)
            elif kind in ('str', 'e', 'b'):
                # preserve non-numeric types so count validation cannot coerce them.
                parsed = value if kind != 'b' else bool(int(value))
            else:
                raise ValueError(f'unsupported cell type: {kind}')
            cells[address] = parsed
        result[name] = cells
    return result


def observations(sheets: dict, selectors: list[dict]) -> list[dict]:
    result = []
    used = set()
    fields = {'construct','condition','sheet','column','batch_column','specimen_column'}
    for select in selectors:
        if set(select) != fields:
            raise ValueError('missing or unknown selection metadata')
        if select['sheet'] not in sheets:
            raise ValueError(f'missing worksheet: {select["sheet"]}')
        cells = sheets[select['sheet']]
        column = select['column']
        if re.fullmatch('[A-Z]+', column) is None:
            raise ValueError('invalid selected column')
        addresses = sorted((a for a in cells if re.fullmatch(column+r'[0-9]+', a)),
                           key=lambda a: int(a[len(column):]))
        selected = 0
        for address in addresses:
            row = int(address[len(column):])
            value = cells[address]
            if row == 1 or value is None:
                continue
            if type(value) is not int or value < 0:
                raise ValueError(f'invalid molecule count at {select["sheet"]}!{address}')
            key = (select['sheet'], address)
            if key in used:
                raise ValueError('observation selected twice')
            used.add(key)
            selected += 1
            result.append({'construct':select['construct'], 'condition':select['condition'],
                           'count':value, 'sheet':select['sheet'], 'cell':address,
                           'batch':cells.get(f'{select["batch_column"]}{row}') if select['batch_column'] else None,
                           'specimen':cells.get(f'{select["specimen_column"]}{row}') if select['specimen_column'] else None,
                           'source':'figure6f.xlsx'})
        if not selected:
            raise ValueError('selected group contains no molecule counts')
    return result


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans('ACGT','TGCA'))[::-1]


def fasta(raw: bytes) -> dict:
    lines = raw.decode('ascii').splitlines()
    match = re.fullmatch(r'>(NC_\d+\.\d+):(\d+)-(\d+) .+', lines[0])
    if not match:
        raise ValueError('expected an accession-version reference interval')
    accession, start, end = match.groups()
    sequence = ''.join(lines[1:])
    if set(sequence)-set('ACGT') or len(sequence) != int(end)-int(start)+1:
        raise ValueError('invalid reference sequence')
    return {'accession':accession, 'start':int(start), 'end':int(end), 'sequence':sequence}


def donors(raw: bytes) -> dict:
    root = xml(office(raw)['word/document.xml'])
    result = {}
    for row in root.iter(W+'tr'):
        cells = [''.join(t.text or '' for t in cell.iter(W+'t')) for cell in row.findall(W+'tc')]
        if len(cells) >= 3 and cells[0].isdigit():
            donor, allele, sequence = cells[:3]
            if donor in result:
                raise ValueError('duplicate donor identifier')
            sequence = sequence.strip().upper()
            if not sequence or set(sequence)-set('ACGT'):
                raise ValueError('invalid donor DNA')
            result[donor] = {'id':donor, 'allele':allele, 'sequence':sequence}
    return result


def matches(reference: dict, query: str) -> list[tuple[int,int]]:
    return [(strand, i) for strand in (1,-1)
            for i in range(len(reference['sequence'])-len(query)+1)
            if (reference['sequence'] if strand == 1 else reverse_complement(reference['sequence'])).startswith(query,i)]


def interval(reference: dict, strand: int, start: int, end: int) -> dict:
    lo, hi = ((reference['start']+start,reference['start']+end-1) if strand==1
              else (reference['end']-end+1,reference['end']-start))
    return {'accession':reference['accession'], 'start':lo, 'end':hi, 'strand':strand}


def reconstruct(donor: dict, reference: dict, contract: dict) -> dict:
    seq = donor['sequence']
    width = contract['homology_arm_bp']
    if len(seq) <= 2*width:
        raise ValueError('donor must contain an insert between its arms')
    left, right = matches(reference,seq[:width]), matches(reference,seq[-width:])
    if len(left) != 1 or len(right) != 1:
        raise ValueError('donor arms must match reference uniquely')
    strand, start = left[0]
    other, end = right[0]
    if strand != other or end < start+width:
        raise ValueError('donor arms have incompatible order or orientation')
    oriented = reference['sequence'] if strand==1 else reverse_complement(reference['sequence'])
    insert = seq[width:-width]
    core = contract['core']
    sites = [(s,i) for s in (1,-1) for i in range(len(insert)-len(core)+1)
             if insert.startswith(core if s==1 else reverse_complement(core),i)]
    if len(sites) != 1:
        raise ValueError('expected exactly one transferred core')
    core_strand, core_start = sites[0]
    return {'donor':donor['id'], 'allele':donor['allele'], 'donor_sequence':seq,
            'reference':oriented[start+width:end], 'insert':insert,
            'reference_bp':end-start-width, 'insert_bp':len(insert),
            'core':core, 'core_start':core_start,
            'core_relative_to_gene':'same' if strand*core_strand==contract['target_gene_strand'] else 'opposite',
            'donor_genomic_strand':strand, 'homology_arm_bp':width,
            'reference_interval':interval(reference,strand,start+width,end),
            'status':'Donor-based intended reconstruction; final clone sequence unverified.'}



def representation(construct: dict, kind: str):
    if kind == 'normalized_core':
        return min(construct['core'], reverse_complement(construct['core']))
    if kind == 'oriented_core':
        return (construct['core'] if construct['core_relative_to_gene']=='same'
                else reverse_complement(construct['core']))
    if kind == 'full_design':
        return [construct.get('donor_sequence'), construct.get('reference'),
                construct['reference_interval']]
    raise ValueError('unknown sequence representation')

def analyze(base: Path = HERE) -> dict:
    case, sources = load_case(base)
    rows = observations(workbook(sources['figure6f.xlsx']), case['groups'])
    for row in rows:
        row['source_sha256'] = sha(sources[row['source']])
    groups = []
    for selection in case['groups']:
        values = [r['count'] for r in rows if (r['construct'],r['condition']) == (selection['construct'],selection['condition'])]
        groups.append({**selection, 'n':len(values), 'sum':sum(values), 'mean':sum(values)/len(values)})
    donor_table = donors(sources['donors.docx'])
    reference, source_reference = fasta(sources['gcy-22.fasta']), fasta(sources['che-1.fasta'])
    constructs = []
    for declared in case['constructs']:
        if declared['donor']:
            item = reconstruct(donor_table[declared['donor']],reference,case['sequence'])
            if len(item['insert']) > len(item['core']):
                found = matches(source_reference,item['insert'])
                if len(found) != 1:
                    raise ValueError('transferred fragment must match its source uniquely')
                strand, start = found[0]
                item['source_interval'] = interval(source_reference,strand,start,start+len(item['insert']))
        else:
            core = case['sequence']['native_core']
            found = matches(reference,core)
            if len(found) != 1:
                raise ValueError('native core must match reference uniquely')
            strand, start = found[0]
            item = {'core':core,'core_relative_to_gene':'same' if strand==case['sequence']['target_gene_strand'] else 'opposite',
                    'reference_interval':interval(reference,strand,start,start+len(core)),
                    'status':'Native gcy-22 promoter in the shared CHE-1 depletion background.'}
        means = {g['condition']:g['mean'] for g in groups if g['construct']==declared['id']}
        item.update(declared)
        item['retention'] = means['auxin']/means['control'] if means['control'] else None
        constructs.append(item)
    comparisons = {}
    for first in constructs:
        for second in constructs:
            if first['id'] == second['id']:
                continue
            def mean(c):
                return next(g['mean'] for g in groups if g['construct']==c['id'] and g['condition']=='auxin')
            delta = mean(second)-mean(first)
            constraints = {}
            for kind in ('normalized_core','oriented_core','full_design'):
                equal = representation(first,kind)==representation(second,kind)
                constraints[kind] = {'equal_inputs':equal, 'minimum_mean_rmse':abs(delta)/2 if equal else 0}
            comparisons[first['id']+':'+second['id']] = {'auxin_mean_difference':delta,'constraints':constraints}
    return {'format':'epigenesis.neural-evidence-report','version':1,'case':case,
            'case_sha256':CASE_SHA256,'implementation_sha256':IMPLEMENTATION_SHA256,
            'observations':rows,'groups':groups,'constructs':constructs,'comparisons':comparisons}


def artifacts(base: Path = HERE) -> dict[str,bytes]:
    report = analyze(base)
    _, sources = load_case(base)
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream,fieldnames=list(report['observations'][0]))
    writer.writeheader()
    writer.writerows(report['observations'])
    encoded = json.dumps(report,ensure_ascii=False,allow_nan=False).replace('<','\\u003c').replace('&','\\u0026')
    template = read_regular_file(HERE/'viewer.html').decode()
    if template.count('__EVIDENCE_DATA__') != 1:
        raise ValueError('viewer must contain exactly one data placeholder')
    output = {'report.json':json_bytes(report), 'observations.csv':stream.getvalue().encode(),
              'index.html':template.replace('__EVIDENCE_DATA__',encoded).encode(),
              'sources/case.json':read_regular_file(base/'case.json')}
    output.update({'sources/'+name:raw for name,raw in sources.items()})
    return output


def write_report(output: Path) -> None:
    files = artifacts()
    output.mkdir(parents=True,exist_ok=False)
    (output/'sources').mkdir()
    for name, raw in files.items():
        with (output/name).open('xb') as target:
            target.write(raw)


def verify_report(output: Path) -> None:
    for name, raw in artifacts().items():
        if read_regular_file(output/name) != raw:
            raise ValueError(f'output differs from source replay: {name}')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('demo','verify'))
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    try:
        if args.command == 'demo':
            write_report(args.output)
        verify_report(args.output)
    except (OSError,ValueError,KeyError,zipfile.BadZipFile,ET.ParseError) as error:
        print(f'evidence: {error}',file=sys.stderr)
        return 1
    print(f'NEURAL-EVIDENCE-VERIFIED: {args.output / "index.html"}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
