# Usage: bash 12.centroMiner.sh namaList

# -r 100 是为了确保鉴定大片段TR的时候程序不堵塞，-l 10000是为了保证小TR不被忽略，后面可以自己筛选
quartet.py CentroMiner -i genome.re.fa --gene genome.renamed.gff -t 30 -r 100 -l 10000

echo -e "Chrid\tStart\tEnd\tLength\tTRlength\tTRcoverage" > genome.centromere.info
length_file="genome.telo.info"
declare -A chr_lengths
while read -r line; do
  chr=$(echo "$line" | cut -f1)
  length=$(echo "$line" | cut -f2)
  chr_lengths[$chr]=$length
done < <(grep -v "#" "$length_file")
chr_length_str=$(for chr in "${!chr_lengths[@]}"; do echo "$chr:${chr_lengths[$chr]}"; done | tr '\n' ',')

centromere_file="./Candidates/*.candidate"
cat $centromere_file | grep -v -e "#" -e "@" | sort -k1,1V -k5,5nr | \
awk -vOFS="\t" -v lengths="$chr_length_str" 'BEGIN {
  split(lengths, chr_lengths, ",");
  for (i in chr_lengths) {
    split(chr_lengths[i], arr, ":");
    length_map[arr[1]] = arr[2];
  }
} {
  if ($2 > 100000 && $3 < length_map[$1] - 100000) {
    print $0;
  }
}' |  awk '!a[$1]++' | awk -vOFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> genome.centromere.info

