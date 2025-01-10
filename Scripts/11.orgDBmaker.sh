# Usage: bash 11.orgDBmaker.sh namaList

cp /public/workspace/biobigdata/project/Plant2t/software/script/11.ko00001.filter4Plant.txt ./

if [ ! -s "11.ko00001.filter4Plant.txt" ]; then
    echo "11.ko00001.filter4Plant.txt not found or empty: run 11.ko00001.filter.R firstly"
    exit 1
fi

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
if [ ! -s "${i}.interpro.gff3" ]; then
    echo "File not found or empty: ${i}.interpro.gff"
    exit 1
fi

if [ ! -s "${i}.pep2ko" ]; then
    echo "File not found or empty: ${i}.pep2ko"
    exit 1
fi

tax_id=$(cat metadata.txt | grep "^Latin Name: " | sed 's/Latin Name: //' | awk -F "[ _\t]" -vOFS=" " '{print $1,$2}' | taxonkit name2taxid | cut -f2)
if [ -z "${tax_id}" ]; then
    tax_id=2759
fi

j=$(cat ${namaList} | cut -f1 | sed 's/[_-]/./g')
if [ ! -f "org.${j}.eg.db_1.0.tar.gz" ]; then
	python3 /public/workspace/biobigdata/project/Plant2t/software/script/11.getGeneTable.py ${i}
	Rscript /public/workspace/biobigdata/project/Plant2t/software/script/11.OrgDBmaker.R ${i} ${j} ${tax_id}
	cp /public/workspace/biobigdata/project/Plant2t/software/script/11.ko00001.filter4Plant.txt ./org.${j}.eg.db/inst/extdata/ko00001.PlanT2T.txt
	bash /public/workspace/biobigdata/project/Plant2t/software/script/11.addPathway2R.sh org.${j}.eg.db
	R CMD build --md5 org.${j}.eg.db
	rm -rf org.${j}.eg.db
fi

done
