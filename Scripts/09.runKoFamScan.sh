# Usage: bash 09.runKoFamScan.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do

cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' > ${i}.pep

exec_annotation -o ${i}.pep2ko.txt \
-p /kofamscan/profiles \
-k kofamscan/ko_list \
--cpu 20 -E 1e-5 -f mapper --report-unannotated ${i}.pep --tmp-dir ./kotmp
sed -e 's/\.[0-9]*//g' ${i}.pep2ko.txt | awk '!a[$1]++' > ${i}.pep2ko

rm ${i}.pep

if [ ! -s "${i}.pep2ko" ]; then
    echo "File empty: ${i}.pep2ko"
    exit 1
else
    echo "Success: ${i}.pep2ko"
fi
    
done    
