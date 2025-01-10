# Usage: bash 16.pepStatic.sh namaList
gffread -y tmp.pep -g genome.re.fa genome.renamed.gff
python3 16.pepstatic.py tmp.pep | cut -f3-5 > pep.static.txt
echo -e "EndStop\tMidStop\tNoStop"
cat pep.static.txt
rm tmp.pep
