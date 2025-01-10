# Usage: bash 20.SaveResultToMySQL.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
python 20.getAnotateSeq.py ${i}
python 20.getGeneFunction.py ${i}
python 20.getGeneTable.py ${i}
bash 20.spe_gene_sql.sh ${i}
done
