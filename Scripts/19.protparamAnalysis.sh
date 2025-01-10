# Usage: bash 19.protparamAnalysis.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do

cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' | sed '/^[^>]/ s/X//g' > ${i}.pep
python /public/workspace/biobigdata/project/Plant2t/software/script/19.protparamAnalysis.py ${i}.pep ${i}.protparam.txt

rm ${i}.pep

done
