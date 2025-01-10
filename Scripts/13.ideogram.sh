# Usage: bash 13.ideogram.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do

j=$(echo ${i} | tr [:upper:] [:lower:] | tr -d '_')

bash 13.ideogramInfo.sh | sed ':a;N;$!ba;s/,\n]/\n]/g' > ${j}.json
if [ ! -s "${j}.json" ]; then
	echo "File empty: ${j}.json"
	exit 1
fi
bash 13.ideogramAnno.sh | sed ':a;N;$!ba;s/\n,/,/g' > ${i}.annotation.json
if [ ! -s "${i}.annotation.json" ]; then
	echo "File empty: ${i}.annotation.json"
	exit 1
fi
cat ${j}.json ${i}.annotation.json

done
