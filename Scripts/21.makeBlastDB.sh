#!/usr/bin/bash
namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read species;do

gffread genome.renamed.gff -g genome.re.fa -x ${species}.cds.fa  -y ${species}.pep.fa

mkdir blast_nDB
mkdir blast_pDB

cut -d " " -f1 genome.re.cds > ./blast_nDB/${species}.cds.fa
cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' > ./blast_pDB/${species}.pep.fa

makeblastdb -in ./blast_nDB/${species}.cds.fa -dbtype nucl
makeblastdb -in ./blast_pDB/${species}.pep.fa -dbtype prot

done
