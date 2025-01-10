# Usage: bash 00.renameGenome.sh namaList

# 从 metadata.txt 提取 Genome File 和 GFF3 File 路径

genome_file=$(cat metadata.txt | grep "Genome File:" | awk -F': ' '{print $2}')
gff3_file=$(cat metadata.txt | grep "GFF3 File:" | awk -F': ' '{print $2}')

# 解压文件
decompress_file() {
    local file=$1
    case $file in
        *.gz)
            gunzip -c "$file" > $(basename "$file" .gz) # 解压并输出为去掉 .gz 后缀的文件名
            echo "$(basename "$file" .gz)"
            ;;
        *.zip)
            unzip -q "$file" -d temp_dir
            real_file=$(ls temp_dir | head -n 1)  # 获取解压后的第一个文件
            mv "temp_dir/$real_file" "$real_file"
            echo "$real_file"
            ;;
        *.tar)
            tar -xf "$file" -C temp_dir
            real_file=$(ls temp_dir | head -n 1)  # 获取解压后的第一个文件
            mv "temp_dir/$real_file" "$real_file"
            echo "$real_file"
            ;;
        *.tar.gz)
            tar -xzf "$file" -C temp_dir
            real_file=$(ls temp_dir | head -n 1)  # 获取解压后的第一个文件
            mv "temp_dir/$real_file" "$real_file"
            echo "$real_file"
            ;;
        *.bz2)
            bunzip2 -c "$file" > $(basename "$file" .bz2) # 解压并输出为去掉 .bz2 后缀的文件名
            echo "$(basename "$file" .bz2)"
            ;;
        *)
            echo "$file" # 非压缩文件直接返回
            ;;
    esac
}

# 解压 genome 文件
genome_real_file=$(decompress_file "$genome_file")
mv "$genome_real_file" genome.fa


# 解压 gff3 文件
gff3_real_file=$(decompress_file "$gff3_file")
mv "$gff3_real_file" genome.gff

# ===================================
# 检查 genome.gff 和 genome.fa 是否匹配, 以及gff3文件是否符合规范
# 1.必须要有gene的feature，2.必须要有ID和Parent属性，3.必须明确正反链，
# 4.feature的起始位置必须小于终止位置，5.gff3中长度必须小于等于fasta文件中的长度，6. 序列名必须一致
python3 00.validate_genome.py genome.gff genome.fa
if [ $? -ne 0 ]; then
    exit 1
fi

# ===================================
# gap number
assembly-stats genome.fa > genome.stats && cat genome.stats
gapNm=$(cat genome.stats | grep -E "^Gaps" | awk -F '[=, ]' '{print $4}')

n90=$(cat genome.stats | grep -E "^N90" | awk -F '[=, ]' '{print $4}')
if [[ $n90 -gt 9000000 ]]; then
    minLength=9000000
else
    minLength=1000000
fi
echo "Selected minLength: $minLength"

faSize -detailed genome.fa > genome.size && cat genome.size

chrNm=$(awk -v minLength=${minLength} '{if ($2 >= minLength){print $1}}' genome.size | wc -l)

# 如果gapNm大于chrNm*2，或者大于50个，就认为这个基因组不符合PlanT2T的要求，退出
if [[ $gapNm -gt $((chrNm*2)) ]] || [[ $gapNm -gt 50 ]]; then
    echo "Too many gaps in the genome, PlanT2T requires that the number of gaps is less than 2 times the number of chromosomes or less than 50."
    exit 1
fi
# ===================================
# 检查蛋白质文件
gffread -y check.pep -g genome.fa genome.gff
if [ ! -f "check.pep" ] || [ ! -s "check.pep" ]; then
    echo "Cannot generate protein sequences from the GFF3 file, please check your GFF3 file."
    exit 1
fi
grep ">" check.pep | cut -d " " -f1 | awk '{ if (a[$1]++) { print "Found duplicate:", $1; exit 1 } }'
if [ $? -ne 0 ]; then
    echo "Exiting due to duplicate sequence names when generating protein sequences"
    exit 1
fi
MidStopAndNoStopFreq=$(python3 16.pepstatic.py check.pep | awk '{print ($4+$5)/$2*100}')
if (( $(echo "$MidStopAndNoStopFreq > 10" | bc -l) )); then
    echo "Too many (${MidStopAndNoStopFreq}% > 10%) sequences with internal stop codons or no stop codons in the protein sequences, please check your GFF3 or Genome file."
    exit 1
fi

# ===================================
# 备份用户提供的原始文件
cp genome.fa genome.raw.fa
cp genome.gff genome.raw.gff

# ===================================
BeforeFilter=$(wc -l genome.size | cut -d " " -f1)
echo "Chromosome number before filtering: $BeforeFilter"

# Filter for sequences that are at least minLength long and save their IDs
awk -v minLength=${minLength} '{if ($2 >= minLength){print $1}}' genome.size | sort -k1,1V > ids.txt && cat ids.txt

# Check the number of sequences after filtering
AfterFilter=$(wc -l ids.txt | cut -d " " -f1)
if [ $AfterFilter -eq 0 ]; then
    echo "Error: All sequences are shorter than the minimum length of $minLength bp."
    exit 1
fi
echo "Chromosome number after filtering: $AfterFilter"

# Extract sequences with the filtered IDs from the genome file and save to a new file
seqkit grep -f ids.txt genome.fa > genome.fa.new

# Filter the GFF3 file to include only entries corresponding to the filtered sequences
awk 'BEGIN{OFS="\t"; FS="\t"} NR==FNR{a[$1]; next} $1 in a {print $0}' ids.txt genome.gff > genome.gff.new

# ===================================
# 查看染色体名字是否以chr开头（不区分大小写），如果不是，则重命名为Chr1,Chr2...
ChrCount=$(grep -c -E "^[Cc][Hh][Rr]" ids.txt)
if [ $ChrCount -ne $AfterFilter ]; then
    echo "Chromosome names do not start with 'Chr', renaming them to 'Chr1', 'Chr2', etc."
    awk -vFS="\t" -vOFS="\t" '{print $1,"Chr"NR}' ids.txt > ids.txt.new
    mv ids.txt.new ids.txt
    # Rename the sequences in the genome file
    awk -vFS="\t" -vOFS="\t" 'NR==FNR{a[$1]=$2; next} /^>/{print ">"a[substr($1,2)],$2; next} {print}' ids.txt genome.fa.new > genome.fa
    # Rename the sequences in the GFF3 file
    awk -vFS="\t" -vOFS="\t" 'NR==FNR{a[$1]=$2; next} $1 in a {$1=a[$1]; print}' ids.txt genome.gff.new > genome.gff
    rm genome.fa.new genome.gff.new
else
    echo "Chromosome names start with 'Chr', no need to rename."
    mv genome.fa.new genome.fa
    mv genome.gff.new genome.gff
fi
# ===================================

faSize -detailed genome.fa > genome.size && cat genome.size

awk -vFS="\t" -vOFS="\t" 'BEGIN{while(getline<"genome.size"){s+=$2}}{print $1,$1,$2,$2/s,"chromosome"}' genome.size > genome.list
rm genome.size ids.txt* check.pep
