#!/bin/bash

# 判断genome.telo.info是否存在
if [ ! -f genome.telo.info ]; then
    echo "genome.telo.info not found, run 3.teloExplorer.sh firstly" >&2
    exit 1
fi

grep -v "#" genome.telo.info | awk -vFS="\t" -vOFS="\t" '{print NR"\t"$0}' | sort -k2,2V | awk -vFS="\t" -vOFS="\t" '
{
    chr_start = 1
    chr_end = $3 - 1
    
    if ($4 == "both") {
        e1 = 1 + $5 * 6
        s2 = $3 - $7 * 6 
        print "Telomere", $2, "#2ca25f", chr_start, e1
        print "Telomere", $2, "#2ca25f", s2, chr_end
    }
    else if ($4 == "right") {
        s1 = chr_end - $7 * 6
        print "Telomere", $2, "#2ca25f", s1, chr_end
    }
    else if ($4 == "left") {
        e1 = chr_start + $5 * 6
        print "Telomere", $2, "#2ca25f", chr_start, e1
    }
}' > Telo.annotation.tmp

grep -v "#" genome.telo.info | awk -vFS="\t" -vOFS="\t" '{print NR"\t"$1}' | sort -k2,2V > genome.telo.info.tmp
awk -vFS="\t" -vOFS="\t" 'NR==FNR {map[$1] = $2; next} $3 in map {print $1, map[$3], "#007bff",$4, $5}' genome.telo.info.tmp ./tmp/genome.telo.label.txt | grep "gap" | sed 's/gap/Gap/g' >> Telo.annotation.tmp

# 初始化JSON结构
echo '{"keys":["name","start","length","color"], "annots":['

# 处理每个染色体的注释
declare -A chr_annots
while IFS=$'\t' read -r type chr color start end; do
    # 跳过注释行
    if [[ $type == "Telomere" || $type == "Gap" ]]; then
        # 计算长度
        length=$((end - start + 1))

        # 添加注释
        annot=$(printf '["%s",%d,%d,"%s"]' "$type" "$start" "$length" "$color")
        chr_annots["$chr"]+="$annot,"
    fi
done < Telo.annotation.tmp

# 合并并输出染色体注释
first=1
for chr in "${!chr_annots[@]}"; do
    if [[ $first -eq 0 ]]; then
        echo ","
    fi
    first=0
    annots="${chr_annots[$chr]}"
    echo "    {\"chr\":\"$chr\",\"annots\":[${annots%,}]}"
done

# 结束JSON结构
echo ']}'



rm genome.telo.info.tmp Telo.annotation.tmp
