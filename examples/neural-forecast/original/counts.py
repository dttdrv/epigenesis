"""Source-bound paired-barcode count contrasts for the declared culture study."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ADMISSION = HERE.parent/'kreimer-replicate-admission-20261018'
sys.path.insert(0,str(ROOT/'examples/neural-regulation'))
import regulation as source

EXACT_INTEGER_LIMIT = 2**(np.finfo(float).nmant+1)


def numeric_row(text, size):
    if re.fullmatch(r'[0-9]+(?:\t[0-9]+)*\r?\n?',text) is None:
        raise ValueError('decimal integer TSV fields required')
    values=np.fromstring(text,sep='\t')
    if values.shape!=(size,):raise ValueError('selected numeric row shape differs')
    return values


def annotation(path, kind):
    with gzip.open(path,'rt') as stream:
        reader=csv.reader(stream,delimiter='\t')
        if next(reader)!=['time','type','rep','bc']:
            raise ValueError('annotation schema differs')
        result={}
        for name,time,type_,rep,barcode in reader:
            key=(int(rep.removeprefix('rep')),int(time.removesuffix('h')),int(barcode))
            if type_!=kind or name in result or key[2]<=0:
                raise ValueError('invalid annotation')
            result[name]=key
    if len(set(result.values()))!=len(result):raise ValueError('duplicate observation coordinate')
    return result


def reduce_pairs(dna, rna, groups, samples):
    if dna.shape!=rna.shape or dna.shape!=groups.shape or dna.ndim!=1 or np.any(groups<0) or np.any(groups>=samples):
        raise ValueError('count/sample shape differs')
    for counts in (dna,rna):
        if not np.isfinite(counts).all() or np.any(counts<0) or np.any(counts>=EXACT_INTEGER_LIMIT) or np.any(counts!=np.floor(counts)):
            raise ValueError('nonnegative exact integer counts required')
    paired=(dna>0)&(rna>0)
    ds=np.bincount(groups[paired],weights=dna[paired],minlength=samples)
    rs=np.bincount(groups[paired],weights=rna[paired],minlength=samples)
    supported=np.bincount(groups[paired],minlength=samples)
    if np.any(ds>=EXACT_INTEGER_LIMIT) or np.any(rs>=EXACT_INTEGER_LIMIT):
        raise ValueError('aggregate outside exact integer range')
    return np.stack((ds,rs,supported))


def contrast(variant, parent):
    variant,parent=np.asarray(variant,dtype=float),np.asarray(parent,dtype=float)
    valid=(variant[...,0]>0)&(variant[...,1]>0)&(parent[...,0]>0)&(parent[...,1]>0)
    result=np.full(valid.shape,np.nan)
    result[valid]=(np.log2(variant[...,1][valid])-np.log2(variant[...,0][valid])
                   -np.log2(parent[...,1][valid])+np.log2(parent[...,0][valid]))
    return result


def missing_samples(rows, variant, parent, indices, cultures, times, admission):
    absent=set(admission['absent_ids'])
    result=[]
    for i in indices:
        for c in cultures:
            for t in times:
                reasons={}
                if variant[i,c,t,2]==0:
                    reasons['variant']='absent_row' if rows[i].id in absent else 'no_matched_positive_slots'
                if parent[i,c,t,2]==0:
                    ids=admission['parent_ids_by_locus'][rows[i].locus]
                    reasons['parent']='all_parent_rows_absent' if all(name in absent for name in ids) else 'no_matched_positive_slots'
                if reasons:result.append({'row':int(i),'culture':c+1,'hour':source.TIMES[t],**reasons})
    return result


def aggregate():
    pins=json.loads((ADMISSION/'DOWNLOAD.json').read_text())
    for item in pins['files']:
        with (ADMISSION/item['name']).open('rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=item['sha256']:
                raise ValueError('raw source identity differs')
    rows,qc=source.load_data()
    rows=[r for r in rows if r.split=='train' and r.kind=='targeted']
    loci={r.locus for r in rows}
    parents={locus:[] for locus in loci}
    for line in (ROOT/'examples/neural-regulation/data/design.fasta').read_text().splitlines():
        if not line.startswith('>WT_region'):continue
        name=line[1:]
        match=re.fullmatch(r'WT_region\d+_(chr[^_:]+:\d+-\d+)',name)
        if match and match[1] in parents and name not in qc['ambiguous_ids']:
            parents[match[1]].append(name)
    parents={k:sorted(set(v)) for k,v in parents.items()}
    if any(not v for v in parents.values()):raise ValueError('no source-qualified WT parent')
    wanted={r.id for r in rows}|{name for names in parents.values() for name in names}
    da=annotation(ADMISSION/'GSE188264_dnaAnnot.tsv.gz','DNA')
    ra=annotation(ADMISSION/'GSE188264_rnaAnnot.tsv.gz','RNA')
    if set(da.values())!=set(ra.values()):raise ValueError('DNA/RNA coordinates differ')
    samples=sorted({key[:2] for key in da.values()})
    cultures=sorted({c for c,t in samples})
    if cultures!=[1,2,3] or samples!=[(c,t) for c in cultures for t in source.TIMES]:
        raise ValueError('assay time/culture contract differs')
    sample_index={key:i for i,key in enumerate(samples)}
    with gzip.open(ADMISSION/'GSE188264_dnaCounts.tsv.gz','rt') as ds,gzip.open(ADMISSION/'GSE188264_rnaCounts.tsv.gz','rt') as rs:
        dh=next(csv.reader([next(ds)],delimiter='\t'));rh=next(csv.reader([next(rs)],delimiter='\t'))
        if len(dh)!=len(da) or set(dh)!=set(da) or len(rh)!=len(ra) or set(rh)!=set(ra):
            raise ValueError('header/annotation mismatch')
        rindex={ra[name]:i for i,name in enumerate(rh)}
        order=np.asarray([rindex[da[name]] for name in dh])
        groups=np.asarray([sample_index[da[name][:2]] for name in dh])
        records={}
        import itertools
        for dline,rline in itertools.zip_longest(ds,rs):
            if dline is None or rline is None:raise ValueError('row count mismatch')
            dn,dtext=dline.split('\t',1);rn,rtext=rline.split('\t',1)
            name=next(csv.reader([dn],delimiter='\t'))[0]
            if next(csv.reader([rn],delimiter='\t'))[0]!=name:raise ValueError('source row alignment differs')
            if name not in wanted:continue
            if name in records:raise ValueError('duplicate selected row')
            d=numeric_row(dtext,len(dh));r=numeric_row(rtext,len(rh))
            records[name]=reduce_pairs(d,r[order],groups,len(samples)).T.reshape(len(cultures),len(source.TIMES),3)
    empty=np.zeros((len(cultures),len(source.TIMES),3))
    variant=np.asarray([records.get(r.id,empty) for r in rows])
    parent=np.asarray([sum((records.get(name,empty) for name in parents[r.locus]),start=empty.copy()) for r in rows])
    if np.any(parent>=EXACT_INTEGER_LIMIT):raise ValueError('pooled parent outside exact integer range')
    effects=contrast(variant,parent)
    report={'rows':len(rows),'families':len({r.group for r in rows}),'cultures':cultures,'times':source.TIMES,
            'source_qualified_ids':len(wanted),'present_ids':len(records),'absent_ids':sorted(wanted-records.keys()),
            'parent_ids_by_locus':parents,
            'parent_totals_by_id':{name:records[name].tolist() for name in sorted({name for names in parents.values() for name in names}) if name in records},
            'rows_with_support_by_culture_time':np.isfinite(effects).sum(axis=0).tolist(),
            'estimator':'log2 pooled RNA/DNA contrast on same-sample matched-positive slots; missing support is NaN',
            'source_scope':'existing TRAIN targeted designs and their source-qualified WT parents only'}
    return rows,effects,variant,parent,report
