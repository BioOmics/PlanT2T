#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import re
import sys
import json


# In[ ]:


species = species = sys.argv[1] if len(sys.argv) > 1 else 'oryza_sativa'
# species = 'arabidopsis_thaliana'
gff_path = './genome.renamed.gff'
function_path = './{}.interpro.gff3'.format(species)
tf_path = './genome.re.pep_output/tf_classification.txt'
kegg_path = './{}.pep2ko'.format(species)

with open(gff_path) as f:
    gff_record = [ rec.strip() for rec in f.readlines() if '#' not in rec ]
    
gene_ID, source_ID = re.findall(r"ID=(.+?);", gff_record[0].split('\t')[8])
# ver_suffix = '.'.join(gene_ID.split('.')[-2:])

function_dict = dict()
go_dict = dict()
try:
    with open(function_path) as f:
        function_record = [ rec.strip() for rec in f.readlines() if '#' not in rec ]
        
    for rec in function_record:
        rec = rec.split('\t')
        if len(rec) != 9:
            continue
        # if ver_suffix not in rec[0]:
        #     rec[0] = rec[0]+'.'+ver_suffix
        # gene = '.'.join([rec[0].split('.')[0]]+rec[0].split('.')[-2:])
        gene = rec[0].split('.')[0]
        function = [rec for rec in re.findall(r'signature_desc=(.+?);', rec[8]) if rec != 'consensus disorder prediction']
        go = re.findall(r'Ontology_term=(.+?);', rec[8])
        go = go[0].replace('"','').split(',') if len(go) != 0 else []
        if gene in function_dict:
            function_dict[gene] += function
        else:
            function_dict[gene] = function

        if gene in go_dict:
            go_dict[gene] += go
        else:
            go_dict[gene] = go

    for gene in function_dict:
        function_dict[gene] = list(set(function_dict[gene]))

    for gene in go_dict:
        go_dict[gene] = list(set(go_dict[gene]))
except:
    print('Error parsing function records. Skip...')

tf_dict = dict()    
try:
    with open(tf_path) as f:
        tf_record = [ rec.strip() for rec in f.readlines() if '#' not in rec ]

    dup_occur = 0
    miss_occur = 0
    for rec in tf_record:
        rec = rec.split('\t')
        # if ver_suffix not in rec[0]:
        #     rec[0] = rec[0]+'.'+ver_suffix
        # gene = '.'.join([rec[0].split('.')[0]]+rec[0].split('.')[-2:])
        gene = rec[0].split('.')[0]
        tf_family = rec[1]
        tf_type = rec[2]
        if gene in tf_dict:
            dup_occur += 1
            if tf_dict[gene]['tf_type'] != tf_type:
                miss_occur += 1
        tf_dict[gene] = {'tf_family': tf_family, 'tf_type': tf_type}
    print("Dup record occur {} times, ratio: {}; miss: {}, ratio: {}; {}!".format(dup_occur, dup_occur/len(tf_record), 
                                                                              miss_occur, miss_occur/len(tf_record),
                                                                       "It's OK" if miss_occur/len(tf_record)<0.1 else "Please Check"))
except:
    print('Error parsing tf records. Skip...')


kegg_dict = dict()
try:
    with open(kegg_path) as f:
        kegg_record = [rec.strip() for rec in f.readlines() if '#' not in rec]
        
    for rec in kegg_record:
        rec = rec.split('\t')
        gene = rec[0].split('.')[0]
        kegg = rec[1] if len(rec) > 1 else ''
        kegg_dict[gene] = kegg
except:
    print('Error parsing KEGG records. Skip...')

# In[ ]:


gff_result = list()
for rec in gff_record:
    rec = rec.split('\t')
    if rec[2] == 'gene':
        gene_ID, source_ID = re.findall(r"ID=(.+?);", rec[8])
        location = rec[0]
        start = int(rec[3]) 
        end = int(rec[4])
        strand = rec[6]
        
    elif rec[2] == 'mRNA':
        if 'Parent' not in rec[8]:
            gene_ID, source_ID = re.findall(r"ID=(.+?);", rec[8])
            location = rec[0]
            start = int(rec[3]) 
            end = int(rec[4])
            strand = rec[6]
        else:
            continue
    else:
        continue
    
    function = ';'.join(function_dict.get(gene_ID, []))
    go = ';'.join(go_dict.get(gene_ID, []))
    tf_family = tf_dict.get(gene_ID, {}).get('tf_family', '')
    tf_type = tf_dict.get(gene_ID, {}).get('tf_type', '')
    kegg = kegg_dict.get(gene_ID, '')

    gff_result.append({'gene_ID': gene_ID, 'source_ID': source_ID, 'go': go, 'kegg': kegg})
    
df = pd.DataFrame(gff_result)
df.to_csv('./{}_gene.csv'.format(species), index=False)

