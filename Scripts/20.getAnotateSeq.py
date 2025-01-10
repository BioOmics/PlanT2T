#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!/usr/bin/python3

## This script is designed for parsing sequences from gff and fasta sequence. Final file will be file which can be imported into SQL database.
## As ortho info need to be loaded in details page, add ortho info here.
### Notice that, intro seq need to start+1 and end-1 when it in annotation type.

from Bio import SeqIO
#from Bio.Seq import Seq
import json
import sys
import subprocess
import re


def get_attriVal(key,attr):
    match = re.search("{}=(.*?);".format(key),attr)
    if match:
        return match.group(1)
    else:
        return 0

species = sys.argv[1] if len(sys.argv) > 1 else 'oryza_sativa'

# Read fasta and gff
with open('genome.renamed.gff') as f:
    gff = f.readlines()

chromSeqs = dict()
for seq in SeqIO.parse('genome.re.fa', 'fasta'):
    chromSeqs[seq.id] = seq.seq

pepSeqs = dict()
for seq in SeqIO.parse('genome.re.pep', 'fasta'):
    pepSeqs[seq.id] = seq.seq


# Build gff dict with needed feature
gff_dict = dict()
for gff_record in gff:
    gff_record = gff_record.strip().split("\t")
    if gff_record[2] in ['gene']:
        final_col = gff_record[-1].split(';')
    else:
        final_col = gff_record[-1].split(';')
        final_col[0],final_col[1] = final_col[1],final_col[0]
    # print(final_col)
    ## 9th col of gff
    parent = ''
    ID = ''
    for col in final_col:
        if col.split('=')[0]=='Parent':
            parent = col.split('Parent=')[1]
            # parent = get_attriVal('Parent',col)
            # print(parent)
        if col.split('=')[0]=='ID':
            ## gff Type
            ID = col.split('=')[1]
            if gff_record[2] == 'mRNA':
                gff_dict[ID] = [{
                    'type': gff_record[2], 
                    'pos': '-'.join([gff_record[0],gff_record[3],gff_record[4]]),
                    'strand': gff_record[6]
                }]
            elif gff_record[2] in ['CDS','five_prime_UTR','three_prime_UTR','exon']:
                ### Note that lines following can be change to meet the fact situation--->Now is ex. LOC_Os01g01610.1
                # print(parent)
                gff_dict[parent].append({
                    'type': gff_record[2],
                    'pos': '-'.join([gff_record[0],gff_record[3],gff_record[4]])
                })
            break
    if gff_record[2] == 'mRNA':
        if parent:
            gff_dict[ID][0]['gene_ID'] = parent
        else:
            gff_dict[ID][0]['gene_ID'] = ID



# In[8]:


result = list()
# Do iteration of gff_dict to get sequence with type
for key in gff_dict:
    gd_record = gff_dict[key]
    gene_id = gd_record[0]['gene_ID']
    subfix = '.'+'.'.join(gene_id.split('.')[-2:])
    gene_chr, gene_s, gene_e = gd_record[0]['pos'].split('-')
    
    ## Get features pos to get intros pos
    feature = list()
    for i in range(1, len(gd_record)):
        _, this_s, this_e = gd_record[i]['pos'].split('-')
        feature.append([int(this_s), int(this_e)])
    feature.sort(key=lambda x:[x[0]])
    
    if len(feature) == 0:
        continue
    
    for i in range(len(feature)):
        if i == 0:
            if feature[i][0]-int(gene_s) > 1:
                gd_record.append({
                    'type':'intro',
                    'pos': gene_chr + '-' + gene_s + "-" + str(feature[i][0]-1)
                })
        else:
            if feature[i][0] - feature[i-1][1] > 1:
                gd_record.append({
                    'type':'intro',
                    'pos': gene_chr + '-' + str(feature[i-1][1]+1) + "-" + str(feature[i][0]-1)
                })
    
    if int(gene_e) - feature[len(feature)-1][1] > 1:
        gd_record.append({
            'type':'intro',
            'pos': gene_chr + '-' + str(feature[len(feature)-1][1]+1) + '-' + gene_e
        })
    
    # Sort all features including intros to meet the output sequence order
    ## Get [1:] to exclude gene record
    strand = gd_record[0]['strand']
    if strand == '+':
        gd_record = sorted(gd_record[1:], key=lambda x:int(x['pos'].split('-')[1]))
    else:
        gd_record = sorted(gd_record[1:], key=lambda x:int(x['pos'].split('-')[1]), reverse=True)
    
    # Get sequence of features
    gene_sequence = list()
    pep_sequence = [{
        'class': 'pep',
        'seq': str(pepSeqs.get(key))
        }]
    if pepSeqs.get(key) is None:
        # print(key)
        continue
    else:
        pep_len = len(pepSeqs.get(key))
    
    transcript_len = 0
    cds_len = 0
    for i in range(len(gd_record)):
        record = gd_record[i]
        if i > 1:
            if 'odd' not in gene_sequence[i-2]['class']:
                record_type = record['type'].lower()+'_odd'
            else:
                record_type = record['type'].lower()
        else:
            record_type = record['type'].lower()+'_odd'
        if 'intro' in record_type:
            record_type = 'intro'
        
        record_pos = record['pos']
        chrom, start, end = record_pos.split("-")
        start, end = [int(start), int(end)]

        ## transcript and cds seq length statistic
        if record_type != 'intro':
            transcript_len += end - start + 1
            if record_type in ['cds','cds_odd']:
                cds_len += end - start + 1
        

        if strand == '+':
            gene_sequence.append({
                'class': record_type, 
                'seq': str(chromSeqs[chrom][start-1:end]
                )})
        else:
            gene_sequence.append({
                'class': record_type, 
                'seq': str(chromSeqs[chrom][start-1:end].reverse_complement()
                )})
    
    tmp_result = {
        'transcript_ID': key,
        'gene_ID': gene_id,
        'genomic_len': int(gene_e)-int(gene_s)+1,
        'transcript_len': transcript_len,
        'cds_len': cds_len,
        'pep_len': pep_len,
        'transcript_location': gene_chr + ':' + gene_s + '-' + gene_e, 
        'gene_strand': strand, 
        'gene_sequence': gene_sequence,
        'pep_sequence': pep_sequence,
        # 'ortho': ortho_dict[key] if key in ortho_dict else ''
    }
    
    result.append(tmp_result)


import pandas as pd
df = pd.DataFrame(result)
df.to_csv('{}_transbase.csv'.format(species), index=False)

