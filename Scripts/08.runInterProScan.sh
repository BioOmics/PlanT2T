# Usage: bash 08.runInterProScan.sh namaList

namaList=$(realpath $1)
cat ${namaList} | cut -f1 | while read i;do

cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' > ${i}.pep
interproscan.sh -i ${i}.pep \
-b ${i}.interpro -f GFF3 \
-cpu 30 --iprlookup --goterms -pa -dp
rm -rf ${i}.pep temp

done
