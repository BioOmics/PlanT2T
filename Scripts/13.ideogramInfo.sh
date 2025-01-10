#!/bin/bash

# 输入文件路径
if [ -d "Candidates" ]; then
	centromere_file="./Candidates/*.candidate"
elif [ -d "candidate" ]; then
	centromere_file="./candidate/*.candidate"
else
	centromere_file="fileNotExit"
fi
length_file="genome.telo.info"

# 判断输入文件是否存在
if [ ! "$(ls -A $centromere_file)" ]; then
  echo "Error: Centromere file not found: $centromere_file, run 12.centroMiner.sh firstly" >&2
  exit 1
fi

if [ ! -f "$length_file" ]; then
  echo "Error: Length file not found: $length_file, run 3.teloExplorer.sh firstly" >&2
  exit 1
fi

# 读取染色体长度信息到一个数组中
declare -A chr_lengths
while read -r line; do
  chr=$(echo "$line" | cut -f1)
  length=$(echo "$line" | cut -f2)
  chr_lengths[$chr]=$length
done < <(grep -v "#" "$length_file")

# 输出格式化的JSON
echo '{"chrBands":['

# 准备染色体长度信息字符串，以供 awk 使用
chr_length_str=$(for chr in "${!chr_lengths[@]}"; do echo "$chr:${chr_lengths[$chr]}"; done | tr '\n' ',')

# 处理着丝粒位置并输出数据
while read -r line; do
  chr=$(echo "$line" | awk '{print $1}')
  cen_start=$(echo "$line" | awk '{print $2}')
  cen_end=$(echo "$line" | awk '{print $3}')
  length=${chr_lengths[$chr]}
  
  # 计算染色体各部分的边界
  cen_center=$(( (cen_start + cen_end) / 2 ))
  
  # 确保正确处理染色体的每一部分
  echo "  \"$chr p p1 1 $cen_start 1 $cen_start gneg\","
  echo "  \"$chr p p1 $cen_start $cen_center $cen_start $cen_center gpos75\","
  echo "  \"$chr q q1 $cen_center $cen_end $cen_center $cen_end gpos75\","
  echo "  \"$chr q q1 $cen_end $length $cen_end $length gneg\","
done < <(cat $centromere_file | grep -v -e "#" -e "@" | sort -k1,1V -k5,5nr | \
awk -vOFS="\t" -v lengths="$chr_length_str" 'BEGIN {
  # 解析染色体长度信息
  split(lengths, chr_lengths, ",");
  for (i in chr_lengths) {
    split(chr_lengths[i], arr, ":");
    length_map[arr[1]] = arr[2];
  }
} {
  # 过滤掉起始位置大于1且末尾位置靠近染色体末尾100000以内的候选着丝粒
  if ($2 > 100000 && $3 < length_map[$1] - 100000) {
    print $0;
  }
}' | awk '!a[$1]++' | awk -vOFS="\t" '{print $1,$2,$3}')

# 关闭JSON格式
echo ']}'
