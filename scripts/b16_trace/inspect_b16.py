import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--data', type=Path, required=True)
parser.add_argument('--evidence', type=Path, required=True)
args = parser.parse_args()
DATA = args.data
OUT = args.evidence
OUT.mkdir(exist_ok=True)

def aligned_bases(fields):
    ref, query = int(fields[3]), 0
    bases = {}
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', fields[5]):
        length = int(length)
        if op in 'M=X':
            for i in range(length):
                bases[ref+i] = (fields[9][query+i], query+i, ord(fields[10][query+i])-33)
            ref += length
            query += length
        elif op in 'DN':
            ref += length
        elif op in 'IS':
            query += length
    return bases

records = []
for line in (DATA / 'b16.combined.sam').read_text().splitlines():
    if not line.startswith('@'):
        f = line.split('\t')
        records.append((f, aligned_bases(f)))

summary = []
for chrom, pos, ref, alt, gene in [('chr4',45802539,'G','C','Aldh1b1'),('chr9',82927102,'G','T','Phip'),('chrX',8125624,'C','A','Wdr13')]:
    rows = []
    for f, bases in records:
        if f[2] != chrom or pos not in bases:
            continue
        base, offset, quality = bases[pos]
        rows.append(dict(read=f[0],flag=int(f[1]),chrom=chrom,start=int(f[3]),mapq=int(f[4]),cigar=f[5],base=base,query_offset_0=offset,base_quality=quality,sequence=f[9],local=f[9][max(0,offset-12):offset]+'['+base+']'+f[9][offset+1:offset+13]))
    with (OUT / f'{gene}-reads.tsv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    fragment_alleles = {}
    for row in rows:
        fragment_alleles.setdefault(row['read'],set()).add(row['base'])
    counts = Counter(row['base'] for row in rows)
    fragment_counts = Counter(next(iter(v)) if len(v)==1 else 'conflict' for v in fragment_alleles.values())
    item = dict(gene=gene,chrom=chrom,pos=pos,ref=ref,alt=alt,reads=len(rows),counts=dict(counts),fragments=len(fragment_alleles),fragment_counts=dict(fragment_counts),flags=dict(Counter(str(r['flag']) for r in rows)),quality_min=min(r['base_quality'] for r in rows),example=next((r for r in rows if r['base']==alt),rows[0]))
    if gene == 'Wdr13':
        supports = [(f,b) for f,b in records if f[2]==chrom and pos in b and b[pos][0]==alt]
        consensus = {}
        for f, bases in supports:
            for p,(base,_,_) in bases.items():
                if abs(p-pos)<=150:
                    consensus.setdefault(p,Counter())[base]+=1
        # A contiguous genomic interval around the mutation, with read support.
        left=pos
        while left-1 in consensus: left-=1
        right=pos
        while right+1 in consensus: right+=1
        sequence=''.join(consensus[p].most_common(1)[0][0] for p in range(left,right+1))
        item['rna_window']=dict(start=left,end=right,genomic_orientation=sequence,mutation_offset_0=pos-left,coverage=[sum(consensus[p].values()) for p in range(left,right+1)])
    summary.append(item)

manifest = json.loads((DATA.parent / 'manifest.json').read_text())
hashes = {}
for name in ['b16.vcf','b16.combined.sam','b16.combined.bam']:
    actual = hashlib.sha256((DATA/name).read_bytes()).hexdigest()
    hashes[name]=actual
    assert actual == manifest['files']['b16.f10/'+name]
(OUT/'summary.json').write_text(json.dumps(dict(summary=summary,sha256=hashes),indent=2))
print(json.dumps(summary,indent=2))
