# Usage: bash 20.SaveResultToMySQL.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
python /public/workspace/biobigdata/project/Plant2t/software/script/20.getAnotateSeq.py ${i}
python /public/workspace/biobigdata/project/Plant2t/software/script/20.getGeneFunction.py ${i}
python /public/workspace/biobigdata/project/Plant2t/software/script/20.getGeneTable.py ${i}
bash /public/workspace/biobigdata/project/Plant2t/software/script/20.spe_gene_sql.sh ${i}
done