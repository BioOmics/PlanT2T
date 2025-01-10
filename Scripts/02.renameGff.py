#!/usr/bin/python3
# coding: utf-8
# Usage: python 02.renameGff.py full_name namaList

import re
import os
import sys
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
import difflib

def stander_fasta(species, dtype, name_dict, c_trans_dict): # dtype need to be fa, cds, pep
    os.chdir("/public/workspace/biobigdata/project/Plant2t/UserUpload")# In condition
    if os.path.exists(species+"/genome.{}".format(dtype)):
        sequenceL = []
        bad_line = []
        for Seq in SeqIO.parse(species+"/genome.{}".format(dtype), "fasta"):
            if name_dict:
                if c_trans_dict:
                    if Seq.id not in c_trans_dict:
                        bad_line.append(Seq.id)
                    else:
                        Seq.id = name_dict[c_trans_dict[Seq.id]]
                        sequenceL.append(Seq)
                else:
                    if Seq.id not in name_dict:
                        bad_line.append(Seq.id)
                    else:
                        Seq.id = name_dict[Seq.id]
                        sequenceL.append(Seq)
            else:
                sequenceL.append(Seq)
        SeqIO.write(sequenceL, species+"/genome.re.{}".format(dtype), "fasta")
        if bad_line:
            with open(species+"/trans_error_line.txt", "w+") as f:
                for i in bad_line:
                    f.write(i+"\n")
    else:
        print(species+"/genome.{} does not exist".format(dtype))
    #os.system("cd {}; mv genome.{} genome.{}.bk; mv genome.re.{} genome.{}; cd -;".format(species,dtype,dtype,dtype,dtype))

def get_attriVal(key,attr):
    match = re.search("{}=(.*?);".format(key),attr)
    if match:
        return match.group(1)
    else:
        return 0

def get_keys(dic, value):
    return [k for k,v in dic.items() if v == value]
    
def judge_line(line): #filter lines do not need
    if line[0] == "#":
        return False
    line = line.strip().split("\t")
    if len(line)!=9:
        return False
    if line[2].lower() not in ["gene","mrna","cds","three_prime_utr","five_prime_utr","transcript","exon"]:
        return False
    else:
        return line[2].lower()
        

def cons_dict(source_record,parent_dict): ## Build 3 levels dict: Gene,Trans,Inter
    result_dict={}
    for index in range(len(source_record)):
        line = source_record[index].split("\t")
        attr = line[8]
        ID = get_attriVal("ID", attr)
        Parent = get_attriVal("Parent", attr)
        if parent_dict==6666 or parent_dict==9999: # Gene_level(6666 has no gene and 9999 has gene)
            if parent_dict == 9999:
                if ID and not Parent:
                    result_dict[index] = ID
            else:
                if line[2] == "mRNA":
                    result_dict[index] = ID
                    line[8] = "ID={};".format(ID) #No Parent attr for mRNA geneline
                    source_record[index] = "\t".join(line)
        else:
            if Parent in parent_dict:
                if Parent in result_dict:
                    result_dict[Parent].append(index)
                else:
                    result_dict[Parent] = [index]
    
    return result_dict
    
    
def sort_record(recordIndex_list, gff_record):
    tmp_list = []
    for index in recordIndex_list:
        line = gff_record[index].split("\t")
        chrom,feature,start = line[0],line[2],line[3]
        tmp_list.append([chrom,feature,int(start),index])
    tmp_list.sort(key=lambda x:(x[0],x[1],x[2]))
    return [x[3] for x in tmp_list]

def renameAttr(target_list, parent, target_type, gff_record, Uniname, UniID, genome_list,type):
    count = 1
    last = 0
    gene_name = {}
    gene_id = {}
    trans_name = {}
    trans_id = {}
    for record in target_list:
        if type == None:
            record_line = gff_record[record].split("\t")
            genome_info = genome_list[record_line[0]]
        else :
            # print(gff_record[record].split("\t")[2])
            if gff_record[record].split("\t")[2] == type:
                # print(gff_record[record])
                record_line = gff_record[record].split("\t")
                # print(record_line)
                genome_info = genome_list[record_line[0]]
            else:
                continue
        if target_type == "Top": # Have ID and no Parent
            if genome_info[2] == "chromosome":
                if last != genome_info[1]:
                    count = 1
                else:
                    count += 1
                last = genome_info[1]
                chr_info = '{:0>2}'.format(genome_info[1].replace("Chr","").replace("chr",""))+"G"
                new_Name = Uniname+'_'+chr_info+'{:0>5}'.format(count)+'0'
                #new_ID = new_Name + ".v1"+".{}".format(UniID)
                new_ID = new_Name
            else:
                if last:
                    count = 1
                    last = 0
                new_Name = Uniname+'_'+genome_info[0]+'{:0>6}'.format(count)+'0'
                new_ID = new_Name + ".v1"+".{}".format(UniID)
                count += 1
            source_ID = get_attriVal("ID", record_line[-1])
            source_Name = get_attriVal("Name", record_line[-1])
            if source_Name:
                gene_name[source_Name] = new_Name# Save gene name
            gene_id[source_ID] = new_Name# Save gene ID
            record_line[-1] = ''.join(["ID="+new_ID,";Name="+new_Name,";Source_ID="+source_ID+";"])
        else:
            parent_Name = get_attriVal("Name", gff_record[parent].split("\t")[-1])
            #new_Parent = parent_Name + ".v1"+".{}".format(UniID) # Same as above
            new_Parent = parent_Name
            if target_type == "Medium": # Have ID and Parent
                source_ID = get_attriVal("ID", record_line[-1])
                source_Name = get_attriVal("Name", record_line[-1])
                source_Accession = get_attriVal("Accession", record_line[-1])
                new_Name = parent_Name + '.{}'.format(count)
                #new_ID = new_Name + ".v1"+".{}".format(UniID)
                new_ID = new_Name
                record_line[-1] = ''.join(["ID=" + str(new_ID), ";Parent=" + str(new_Parent), ";Name=" + str(new_Name), ";Source_ID=" + str(source_ID) + ";"])
                if source_Accession:
                    trans_name[source_Accession] = new_Name
                elif source_Name:
                    trans_name[source_Name] = new_Name
                trans_id[source_ID] = new_Name
                    
            elif target_type == "Bottom": # Have Parent an no ID
                new_Name = parent_Name + "." + record_line[2] + '.{}'.format(count)
                #new_ID = new_Name + ".v1"+".{}".format(UniID)
                new_ID = new_Name
                source_ID = get_attriVal("ID", record_line[-1])
                if source_ID:
                    #record_line[-1] = "Parent="+new_Parent+";ID="+new_ID+";Source_ID="+source_ID+";"
                    record_line[-1] = "ID="+ new_ID +";Parent="+new_Parent+";Source_ID=" + source_ID + ";"
                else:
                    record_line[-1] = "ID="+new_ID+";Parent="+new_Parent+";"
            count += 1
        
        gff_record[record] = '\t'.join(record_line)
        
    if target_type == "Top":
        return gene_name, gene_id
    elif target_type == "Medium":
        return trans_name, trans_id
        

def gff_output(file_name, dicts):
    levels = len(dicts)
    if levels == 3:
        sorted_Gene_list, Trans_dict, Inter_dict = dicts
    else:
        sorted_Gene_list, Trans_dict = dicts
    with open(file_name,"w+") as f:
        for ID in sorted_Gene_list:
            gene_line = filtered_gff[ID]
            f.write(gene_line+"\n")
            for sub_ID in Trans_dict[Gene_dict[ID]]:
                mRNA_line = filtered_gff[sub_ID]
                f.write(mRNA_line+"\n")
                if levels == 3:
                    inter = get_attriVal("Source_ID", filtered_gff[sub_ID].split("\t")[-1])
                    if inter in Inter_dict:
                        for ssub_ID in Inter_dict[get_attriVal("Source_ID", filtered_gff[sub_ID].split("\t")[-1])]:
                            inter_line = filtered_gff[ssub_ID]
                            f.write(inter_line+"\n")
        
def check_name_list(name_dict, source_name):
    match = None
    rev_flag = False
    trans_name = {}
    
    
    for i in name_dict:
        first = i
        break
    
    if first not in source_name:
        match = difflib.get_close_matches(first, source_name, n=1, cutoff=0.2)
        if match:
            match = match[0]
            pattern_match = ''.join(i[-1] for i in difflib.ndiff(first, match) if i[0] == "+")
            pattern_self = ''.join(i[-1] for i in difflib.ndiff(match, first) if i[0] == "+")
        else:
            return trans_name

        for i in name_dict:
            if pattern_self:
                parsed = i.replace(pattern_self, pattern_match)
                #if parsed in source_name:
                trans_name[parsed] = i
            else:
                rev_flag = True
                   
    if rev_flag:
        for i in source_name:
            parsed = i.replace(pattern_match, pattern_self)
            #if parsed in name_dict_key:
            trans_name[i] = parsed        
            
    return trans_name


os.chdir("/public/workspace/biobigdata/project/Plant2t/UserUpload/")
species = sys.argv[1]
name_list = sys.argv[2]
#species = 'aegilops_tauschii'

if not species:
    print("Please enter speciesName")
    exit()
    
if not name_list:
    print("Please enter nameList")
    exit()
    
genome_list = open("{}/genome.list".format(species)).readlines()
genome_list_dict = {}
for line in genome_list:
    line=line.strip().split("\t")
    genome_list_dict[line[0]] = [line[0],line[1],line[-1]]
    

name_list = open(name_list).readlines()
name_dict = {}
for i in range(len(name_list)):
    name=name_list[i].strip().split()
    name_dict[name[0]]=[name[1],i]

## Step 1 - Preparse
print("Start standardization of '{}' fasta format".format(species))
stander_fasta(species,"fa",None,None)
print("Fasta standardization done. GFF parsing:")

gff_file = open("{}/genome.gff".format(species)).readlines()
filtered_gff = []
feature_type = set()
for index in range(len(gff_file)):
    line = gff_file[index]
    feature = judge_line(line)
    if feature:
        feature_type.add(feature)
        filtered_gff.append(gff_file[index].strip()+";")# Some attr has no ; in $
        # filtered_gff.append(gff_file[index].strip())
        
if "gene" not in feature_type:
    gene_flag = 6666
else:
    gene_flag = 9999

print("Start construction of dicts")
start_time = time() # Starting time of cons dict
Gene_dict = cons_dict(filtered_gff, gene_flag)#gene所在行：geneID
tGene_dict = {Gene_dict[key]:key for key in Gene_dict}#geneID:gene所在行
Trans_dict = cons_dict(filtered_gff, tGene_dict)#geneID：[mRNA所在行（可能有多个剪切）]
tTrans_dict = {}
tTrans_dict_Gname = {}
for i in Trans_dict:
    for j in Trans_dict[i]:
        trans_ID = get_attriVal("ID",filtered_gff[j].split("\t")[8])
        trans_feature = filtered_gff[j].split("\t")[2]
        if trans_feature != "CDS":
            tTrans_dict[trans_ID] = j#mRNAID:mRNA所在行
            tTrans_dict_Gname[trans_ID] = i#mRNAID:geneID
if len(tTrans_dict) > 1:
    Inter_dict = cons_dict(filtered_gff, tTrans_dict)#mRNAID:[CDS所在行]
else:
    Inter_dict =None
end_time = time() # Endding time of cons dict
if Inter_dict: # Judging levels of dicts
    levels = 3
else:
    levels = 2
print("{} dicts construction done. Levels {}, time {}s.".format(species, levels, end_time-start_time))



## Sort all record in annotation tree
Gene_index = [ key for key in Gene_dict ]#gene所在行
sorted_Gene_list = sort_record(Gene_index, filtered_gff)#按gene的strat从小到大排序所在行

no_value = []
no_CDS = []
for gindex in sorted_Gene_list:
    if Gene_dict[gindex] in Trans_dict:
        Trans_dict[Gene_dict[gindex]] = sort_record(Trans_dict[Gene_dict[gindex]], filtered_gff)
        if levels == 3:
            for tindex in Trans_dict[Gene_dict[gindex]]:
                transID = get_attriVal("ID",filtered_gff[tindex].split("\t")[8])
                if transID in Inter_dict:
                    Inter_dict[transID] = sort_record(Inter_dict[transID], filtered_gff)
                else:
                    no_CDS.append(transID)
    else:
        no_value.append(gindex)

#Check if no_value conclude no_CDS
if no_CDS:
    for i in no_CDS:
        if i in tTrans_dict:
            normal_Genekey = get_keys(Trans_dict,[tTrans_dict[i]])#找geneID
            if normal_Genekey:
                no_value.append(tGene_dict[normal_Genekey[0]])#找gene所在行
            else:
                Trans_dict[tTrans_dict_Gname[i]].pop(Trans_dict[tTrans_dict_Gname[i]].index(tTrans_dict[i]))
                if not Trans_dict[tTrans_dict_Gname[i]]:
                    no_value.append(tGene_dict[tTrans_dict_Gname[i]])

#Del value in no_value
for i in range(len(sorted_Gene_list)):
    if sorted_Gene_list[i] in no_value:
        Gene_dict.pop(sorted_Gene_list[i])
tGene_dict = {Gene_dict[key]:key for key in Gene_dict}
for i in no_value:
    sorted_Gene_list.pop(sorted_Gene_list.index(i))
    

## Step2 - Rename            
## Gene level
gene_name, gene_id = renameAttr(sorted_Gene_list, None, "Top", filtered_gff, name_dict[species][0], name_dict[species][1],
          genome_list_dict,None)#Top:gene level, no parent, have ID
if levels == 3:
    max_Trans_dict = {}
    trans_name = {}
    trans_id = {}
for ID in tGene_dict:
    if levels == 3:
        trans_name_tmp, trans_id_tmp = renameAttr(Trans_dict[ID], tGene_dict[ID], "Medium", filtered_gff, None, name_dict[species][1], genome_list_dict,None)
        for key in trans_name_tmp:
            trans_name[key] = trans_name_tmp[key]
        for key in trans_id_tmp:
            trans_id[key] = trans_id_tmp[key]
        max_trans_length = 0 ## longest trans
        for sub_ID in Trans_dict[ID]:
            trans_line = filtered_gff[sub_ID].split("\t")
            length = int(trans_line[4])-int(trans_line[3])
            if length > max_trans_length:
                max_trans_length = length
                max_Trans_dict[ID] = [sub_ID]
            source_ID = get_attriVal("Source_ID", trans_line[-1])
            if source_ID in Inter_dict:
                # print(Inter_dict[source_ID])
                renameAttr(Inter_dict[source_ID], sub_ID, "Bottom", filtered_gff, None, name_dict[species][1], genome_list_dict, "CDS")
                renameAttr(Inter_dict[source_ID], sub_ID, "Bottom", filtered_gff, None, name_dict[species][1],genome_list_dict, "exon")
                renameAttr(Inter_dict[source_ID], sub_ID, "Bottom", filtered_gff, None, name_dict[species][1],genome_list_dict, "three_prime_UTR")
                renameAttr(Inter_dict[source_ID], sub_ID, "Bottom", filtered_gff, None, name_dict[species][1],genome_list_dict, "five_prime_UTR")
    else:
        renameAttr(Trans_dict[ID], tGene_dict[ID], "Bottom", filtered_gff, None, name_dict[species][1], genome_list_dict,None)


## Merge name dict
name_dict = {}
if levels == 3:
    name_dict.update(trans_name)
    name_dict.update(trans_id)
else:
    name_dict.update(gene_name)
    name_dict.update(gene_id)
    
        
## imputation
if levels == 3:
    for i in Trans_dict:
        if i in max_Trans_dict:
            pass
        else:
            max_Trans_dict[i]=Trans_dict[i]
        

## Step3 - New gff output
os.chdir("/public/workspace/biobigdata/project/Plant2t/UserUpload/{}".format(species))
if levels == 3:
    print("Outputing 2 files...")
    gff_output("genome.renamed.gff", [sorted_Gene_list, Trans_dict, Inter_dict])
    gff_output("genome.longest.gff", [sorted_Gene_list, max_Trans_dict, Inter_dict])
else:
    print("Outputing 1 file...")
    gff_output("genome.renamed.gff", [sorted_Gene_list, Trans_dict])
print("Outputing Done.")
    
## Step4 - Change name from name dict
print("Renaming cds and pep start")
os.chdir("/public/workspace/biobigdata/project/Plant2t/UserUpload/{}".format(species))

source_cds_name = []
source_pep_name = []

for Seq in SeqIO.parse("genome.cds", "fasta"):
    source_cds_name.append(Seq.id)
for Seq in SeqIO.parse("genome.pep", "fasta"):
    source_pep_name.append(Seq.id)  
    

print("Start name_list check of genome.cds and genome.pep...")
start_time = time()
c_name = [None]*2
c_name[0] = check_name_list(name_dict, source_cds_name)
c_name[1] = check_name_list(name_dict, source_pep_name)
# print(c_name[0])
end_time = time()

if not c_name[0]:
    print("Stander CDS name list, no change.")
else:
    print("Change CDS name list, time:{}.".format(end_time-start_time))
if not c_name[1]:
    print("Stander PEP name list, no change.")
else:
    print("Change PEP name list, time:{}.".format(end_time-start_time))

dtype = ["cds","pep"]
for num in range(len(dtype)):
    try:
        stander_fasta(species, dtype[num], name_dict, c_name[num])
    except:
        print("{} {} trans_error".format(species, dtype[num]))
        os.system("echo {} >> /public/workspace/biobigdata/project/Plant2t/UserUpload/trans_error.log".format(species))

with open("{}/all_name.py".format(species),"w+") as f:
    if levels == 3:
        f.write("trans_name=")
        f.write(str(trans_name)+"\n")
    f.write("gene_name=")
    f.write(str(gene_name)+"\n")
         
print("Renaming Done.")
print("Species {} Done.\n".format(species))