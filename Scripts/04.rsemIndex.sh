# Usage: bash 04.rsemIndex.sh namaList

# mkdir RSEMIndex
# cd RSEMIndex

# rsem-gff3-to-gtf ../genome.renamed.gff ./genome.renamed.gtf
# rsem-prepare-reference --gtf ./genome.renamed.gtf ../genome.re.fa ./RSEMIndex -p 30 --star

# if [[ ! -f "./RSEMIndex.n2g.idx.fa" ]];then
#     echo "Fail to build RSEM index"
#     exit 1
# fi

echo "Skip RSEM and STAR index building"

