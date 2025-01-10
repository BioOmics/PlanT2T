#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json
import re
import sys
import os

species = sys.argv[1] if len(sys.argv) > 1 else 'oryza_sativa'
tfend = {}
converted_data = {}
js_data = []
if os.path.exists('{}.interpro.gff3'.format(species)):
    with open('{}.interpro.gff3'.format(species), 'r') as file:
        for line in file:
            if line.startswith('##sequence-region'):
                items= line.strip().split(' ')
                id = items[1]
                end = int(items[3])
                tfend[id] = end
        current_Trans_ID = None
        current_Gene_ID = None
    with open('{}.interpro.gff3'.format(species), 'r') as file:
        for line in file:
            if line.startswith('##') :
                continue
            fields = line.strip().split('\t')
            if len(fields)  < 9 :
                continue
            if len(fields) == 9 and fields[1] != '.' :
                seqid = fields[0]
                source = fields[1]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                key_value_pairs = attributes.split(';')
                attributes_dict = {}
                for pair in key_value_pairs:
                    if '=' in pair:
                        key, value = pair.split('=',1)
                        attributes_dict[key] = value
                    else:
                        continue
                Name = attributes_dict.get('Name')
                signature_desc = attributes_dict.get('signature_desc')
                signature_desc = signature_desc if signature_desc != 'consensus disorder prediction' else None
                Dbxref = re.findall(r'"InterPro:(.+?)"', attributes_dict.get('Dbxref')) if attributes_dict.get('Dbxref') else []
                Dbxref = Dbxref[0] if len(Dbxref) != 0 else None
                Trans_ID = seqid
                # Gene_ID = seqid.split('.')[0]+'.'+'.'.join(seqid.split('.')[-2:])
                Gene_ID = seqid.split('.')[0]
                tfend1 = tfend.get(seqid,0)
                js_obj = {
                    'transcript_ID': Trans_ID,
                    'gene_ID': Gene_ID,
                    'source': source,
                    'start': start,
                    'end': end,
                    'name': Name,
                    'signature_desc': signature_desc,
                    'Dbxref': Dbxref,
                    'tfend': tfend1
                }

                js_data.append(js_obj)
    #             if Gene_ID in converted_data:
    #                 if Trans_ID in converted_data[Gene_ID]:
    #                     converted_data[Gene_ID][Trans_ID].append(js_obj)
    #                 else:
    #                     converted_data[Gene_ID][Trans_ID] = [js_obj]
    #             else:
    #                 converted_data[Gene_ID] = {Trans_ID: [js_obj]}

    # converted_data_json = json.dumps(converted_data, indent=2)

    # with open('./data/{}.geneFunction.json'.format(species), 'w') as json_file:
    #     json_file.write(converted_data_json)
else:
    print("No interpro annotation. Skip...")


# In[ ]:


import pandas as pd
df = pd.DataFrame(js_data)
df.to_csv('{}_transfunc.csv'.format(species))

