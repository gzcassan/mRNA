import argparse
import csv
import gzip
import json
import math
import struct
from collections import Counter
from pathlib import Path
from Bio.Seq import Seq

parser=argparse.ArgumentParser()
parser.add_argument('--data', type=Path, required=True)
parser.add_argument('--evidence', type=Path, required=True)
args=parser.parse_args()
out=args.evidence
data=args.data
# Decode original BAM records and compare with the readable SAM fixture.
# BAM stores each SEQ base in four bits; CIGAR operations are packed integers.
raw=gzip.decompress((data/'b16.combined.bam').read_bytes())
assert raw[:4]==b'BAM\x01'
offset=4
ltext=struct.unpack_from('<i',raw,offset)[0]; offset+=4+ltext
nref=struct.unpack_from('<i',raw,offset)[0]; offset+=4
refs=[]
for _ in range(nref):
    size=struct.unpack_from('<i',raw,offset)[0];offset+=4
    refs.append(raw[offset:offset+size-1].decode());offset+=size+4
decoded=[]
while offset<len(raw):
    size=struct.unpack_from('<i',raw,offset)[0];offset+=4
    block=raw[offset:offset+size];offset+=size
    rid,pos,bin_mq_nl,flag_nc,lseq,nrid,npos,tlen=struct.unpack_from('<iiIIiiii',block)
    lname=bin_mq_nl&255;mapq=(bin_mq_nl>>8)&255
    nc=flag_nc&65535;flag=flag_nc>>16
    name=block[32:32+lname-1].decode()
    cigars=struct.unpack_from('<'+'I'*nc,block,32+lname)
    cigar=''.join(str(c>>4)+'MIDNSHP=XB'[c&15] for c in cigars)
    seqbytes=block[32+lname+4*nc:32+lname+4*nc+(lseq+1)//2]
    seq=''.join('=ACMGRSVTWYHKDBN'[(seqbytes[i//2]>>(4 if i%2==0 else 0))&15] for i in range(lseq))
    decoded.append((name,str(flag),refs[rid],str(pos+1),str(mapq),cigar,seq))
sam=[]
for line in (data/'b16.combined.sam').read_text().splitlines():
    if line.startswith('@'): continue
    f=line.split('\t');sam.append(tuple(f[i] for i in [0,1,2,3,4,5,9]))
assert Counter(decoded)==Counter(sam)
print('BAM-SAM sequence/alignment parity:',len(decoded),'records')

evidence=json.loads((out/'summary.json').read_text())
wdr=evidence['summary'][2]
window=wdr['rna_window']
reference=json.loads((out/'wdr13-mm10-full-gene.json').read_text())
genome=reference['dna'].upper();start=reference['start']
genes=json.loads((out/'wdr13-mm10-genes.json').read_text())['ncbiRefSeq']
gene=next(g for g in genes if g['name']=='NM_026137.5')
positions=[]
for left,right in zip(map(int,gene['exonStarts'].strip(',').split(',')),map(int,gene['exonEnds'].strip(',').split(','))):
    positions.extend(range(max(left,gene['cdsStart']),min(right,gene['cdsEnd'])))
positions.reverse()
cds=''.join(str(Seq(genome[p-start]).complement()) for p in positions)
protein=str(Seq(cds).translate())
mut_index=positions.index(8125624-1)
mut_cds=cds[:mut_index]+str(Seq('A').complement())+cds[mut_index+1:]
mut_protein=str(Seq(mut_cds).translate())
aa_index=mut_index//3
assert protein[aa_index]=='S' and mut_protein[aa_index]=='I'
local_positions=positions[(aa_index-12)*3:(aa_index+13)*3]
rna_cds=''.join(str(Seq(window['genomic_orientation'][p+1-window['start']]).complement()) for p in local_positions)
rna_protein=str(Seq(rna_cds).translate())
assert rna_protein=='KLQGHSAPVLDVIVNCDESLLASSD'
wt_protein=protein[aa_index-12:aa_index+13]
diffs=[]
for p in range(window['start'],window['end']+1):
    readbase=window['genomic_orientation'][p-window['start']]
    refbase=genome[p-1-start]
    if readbase!=refbase: diffs.append(dict(pos=p,ref=refbase,rna=readbase))
codon_start=mut_index//3*3
translation=dict(transcript=gene['name'],strand=gene['strand'],aa_position_1=aa_index+1,cds_position_1=mut_index+1,reference_codon=cds[codon_start:codon_start+3],mutant_codon=mut_cds[codon_start:codon_start+3],reference_protein_window=wt_protein,rna_protein_window=rna_protein,reference_rna_differences=diffs,window_base_support_min=min(window['coverage'][p+1-window['start']] for p in local_positions),bam_sam_records_verified=len(decoded))
(out/'translation.json').write_text(json.dumps(translation,indent=2))
print(json.dumps(translation,indent=2))

rows=list(csv.DictReader((out/'iedb-ba.tsv').open(),delimiter='\t'))
assert len(rows)==396
mutant=[r for r in rows if r['seq_num']=='1' and int(r['start'])<=13<=int(r['end'])]
assert len(mutant)==76
comparison=[]
for r in sorted(mutant,key=lambda r:float(r['ic50'])):
    wt=next(w for w in rows if w['seq_num']=='2' and w['allele']==r['allele'] and w['start']==r['start'] and w['end']==r['end'])
    reference_row=next(w for w in rows if w['seq_num']=='3' and w['allele']==r['allele'] and w['start']==r['start'] and w['end']==r['end'])
    value=float(r['ic50'])
    score=(1+math.exp(-350/150))/(1+math.exp((value-350)/150)) if value<5000 else 0
    comparison.append(dict(allele=r['allele'],start=r['start'],end=r['end'],mutant_peptide=r['peptide'],mutant_ic50_nM=r['ic50'],mutant_BA_rank_percent=r['percentile_rank'],reverted_background_peptide=wt['peptide'],reverted_background_ic50_nM=wt['ic50'],reverted_background_BA_rank_percent=wt['percentile_rank'],mm10_reference_peptide=reference_row['peptide'],mm10_reference_ic50_nM=reference_row['ic50'],mm10_reference_BA_rank_percent=reference_row['percentile_rank'],legacy_affinity_score=score,exploratory_rank_le_2=float(r['percentile_rank'])<=2,strict_ic50_le_500=value<=500))
with (out/'Wdr13-peptide-comparisons.tsv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(comparison[0]),delimiter='\t');writer.writeheader();writer.writerows(comparison)
passing=[r for r in comparison if float(r['mutant_ic50_nM'])<5000]
score=sum(r['legacy_affinity_score'] for r in passing)
ranking=dict(mutant_containing_pairs=len(comparison),legacy_affinity_cutoff_5000_passing=len(passing),exploratory_rank_2_passing=sum(r['exploratory_rank_le_2'] for r in comparison),strict_affinity_500_passing=sum(r['strict_ic50_le_500'] for r in comparison),legacy_target_epitope_score=score,expression_score=math.sqrt(21),illustrative_combined_score=math.sqrt(21)*score,scope='manual single-window recalculation; not a Vaxrank CLI run or final mRNA selection')
(out/'ranking.json').write_text(json.dumps(ranking,indent=2))
print('Top comparisons',json.dumps(comparison[:8],indent=2))
print(json.dumps(ranking,indent=2))
