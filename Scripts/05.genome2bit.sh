# Usage: bash 05.genome2bit.sh namaList

faToTwoBit genome.re.fa genome.re.2bit
python 05.calculateGCcontent.py genome.re.fa > genome.re.GC.txt
cat genome.re.GC.txt
