# Usage: bash 14.getKEGG.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
species=$i
echo ${species}

cp 11.ko00001.filter4Plant.txt ./

awk -vFS="\t" -vOFS="\t" '{print $5, $3}' 11.ko00001.filter4Plant.txt > level3.txt
awk -vFS="\t" -vOFS="\t" '{print $5, $2}' 11.ko00001.filter4Plant.txt > level2.txt
awk -vFS="\t" -vOFS="\t" '{print $5, $1}' 11.ko00001.filter4Plant.txt > level1.txt

awk -vFS="\t" -vOFS="\t" '$2 != "" {print $2, $1}' ./${species}.pep2ko > level.gene2ko.txt

join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 level3.txt) <(sort -t $'\t' -k1,1 level.gene2ko.txt) > level3.merged.txt
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 level2.txt) <(sort -t $'\t' -k1,1 level.gene2ko.txt) > level2.merged.txt
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 level1.txt) <(sort -t $'\t' -k1,1 level.gene2ko.txt) > level1.merged.txt

awk -F'\t' 'BEGIN {OFS="\t"} {print "P", $2, $3}' level3.merged.txt > level3.final.txt
awk -F'\t' 'BEGIN {OFS="\t"} {print "C", $2, $3}' level2.merged.txt > level2.final.txt
awk -F'\t' 'BEGIN {OFS="\t"} {print "F", $2, $3}' level1.merged.txt > level1.final.txt

cat level*.final.txt | awk '!a[$1$2$3]' > kegg_term.txt

python3 14.getKEGGnpy.py
mv kegg_dict.npy ${species}_kegg.npy
rm kegg_term.txt level*.txt 11.ko00001.filter4Plant.txt

done
