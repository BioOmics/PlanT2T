# Usage: bash 02.renameGff.sh namaList
namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
	python /public/workspace/biobigdata/project/Plant2t/software/script/02.renameGff.py ${i} ${namaList}
	if [ ! -f "genome.re.pep" ]; then
		echo "File not found or empty: genome.re.pep"
		exit 1
	fi
	head -2 genome.re.pep
	grep ">" genome.re.pep | cut -d " " -f1 | awk '{ if (a[$1]++) { print "Found duplicate:", $1; exit 1 } }'
	if [ $? -ne 0 ]; then
		echo "Exiting due to duplicate sequence names"
		exit 1
	fi
done
