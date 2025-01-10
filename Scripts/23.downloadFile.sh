#!/usr/bin/bash

namaList=$(realpath $1)

echo "Compressing files"

mkdir download

cat ${namaList} | cut -f1 | while read i;do
    mv genome.re.fa ./download/${i}.genome.fa && gzip ./download/${i}.genome.fa
    mv genome.renamed.gff ./download/${i}.gff3 && gzip ./download/${i}.gff3
    cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' > ./download/${i}.pep.fa && gzip ./download/${i}.pep.fa
    cut -d " " -f1 genome.re.cds > ./download/${i}.cds.fa && gzip ./download/${i}.cds.fa
    mv genome.raw.fa ./download/${i}.raw.fa && gzip ./download/${i}.raw.fa
    mv genome.raw.gff ./download/${i}.raw.gff3 && gzip ./download/${i}.raw.gff3
    cp genome.telo.info ./download/${i}.telomere.txt
    cp genome.centromere.info ./download/${i}.centromere.txt
done