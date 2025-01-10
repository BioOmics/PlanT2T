#!/usr/bin/python3

import numpy as np
with open('kegg_term.txt') as f:
   content = [ rec.strip() for rec in f.readlines() ]

kegg_dict = dict()
for rec in content:
    t,term,gene = rec.split('\t')
    if t not in kegg_dict:
        kegg_dict[t] = dict()
    if term not in kegg_dict[t]:
        kegg_dict[t][term] = list()
    kegg_dict[t][term].append(gene)

np.save('kegg_dict.npy', kegg_dict)
