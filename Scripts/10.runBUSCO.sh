# Usage: bash 10.runBUSCO.sh namaList

thread=20
/public/workspace/biobigdata/project/Plant2t/software/busco-5.8.2/bin/busco -i genome.re.fa -o busco_result -m genome --offline -c ${thread} -f \
-l /public/workspace/biobigdata/project/Plant2t/software/busco-5.8.2/buscodb/embryophyta_odb10