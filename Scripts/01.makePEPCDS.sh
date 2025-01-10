# Usage: bash 01.makePEPCDS.sh namaList

# Error: discarding overlapping duplicate gene feature
# awk '!a[$1$3$4$5]++' genome.gff > genome.dup.gff
# mv genome.gff genome.remove.gff && mv genome.dup.gff genome.gff
gffread -x genome.cds -y genome.pep -g genome.fa genome.gff

# check if the genome.cds and genome.pep exists and is not empty
if [ ! -f "genome.cds" ] || [ ! -s "genome.cds" ]; then
    echo "File not found or empty: genome.cds"
    exit 1
else
    echo "genome.cds is ready."
fi

if [ ! -f "genome.pep" ] || [ ! -s "genome.pep" ]; then
    echo "File not found or empty: genome.pep"
    exit 1
else
    echo "genome.pep is ready."
fi
