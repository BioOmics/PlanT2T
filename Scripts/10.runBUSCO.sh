# Usage: bash 10.runBUSCO.sh namaList

thread=20
busco -i genome.re.fa -o busco_result -m genome --offline -c ${thread} -f \
-l embryophyta_odb10
