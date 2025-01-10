#!/usr/bin/bash
namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read genome;do

wsdir="/public/workspace/plant2tUserGenomes/${genome}"
jbdir="/public/workspace/jbrowser2"

mkdir -p ${wsdir}
cp genome.re.fa ${wsdir}/${genome}.fa
cp genome.renamed.gff ${wsdir}/${genome}.gff3
samtools faidx ${wsdir}/${genome}.fa

config_json='{"formatDetails": {"feature":"jexl: {Plant2t: \"<a style=text-decoration:none;color:#2ca25f;background-color:#e5f5f9;padding:4px;border-radius:5px; href=https://biobigdata.nju.edu.cn/plant2t/genome/\"+feature.dbxref+\"/details?gene=\"+feature.id+\">GO</a>\",dbxref:undefined}", "subfeatures":"jexl: {Plant2t: \"<a style=text-decoration:none;color:#2ca25f;background-color:#e5f5f9;padding:4px;border-radius:5px; href=https://biobigdata.nju.edu.cn/plant2t/genome/\"+feature.dbxref+\"/details?gene=\"+feature.parent+\"&trans=\"+feature.id+\">GO</a>\",dbxref:undefined}"}}'

cd ${jbdir}
npx @jbrowse/cli add-assembly ${wsdir}/${genome}.fa --load symlink --out ./jbrowse2/plant2tupload/${genome}
npx @jbrowse/cli sort-gff ${wsdir}/${genome}.gff3 | awk -vFS="\t" -vOFS="\t" -vdb=${genome} '{if($1!~/#/){print $0"dbxref="db";"}else{print $0}}' | bgzip > ${wsdir}/${genome}.gff3.gz
tabix ${wsdir}/${genome}.gff3.gz
npx @jbrowse/cli add-track ${wsdir}/${genome}.gff3.gz --load symlink --out ./jbrowse2/plant2tupload/${genome} --config "${config_json}"
npx @jbrowse/cli text-index --out ./jbrowse2/plant2tupload/${genome} --tracks=${genome}.gff3

if [ $(zcat ${wsdir}/${genome}.gff3.gz | awk '($3=="gene"){print $1}' | wc -l) -gt 0 ];then
    zcat ${wsdir}/${genome}.gff3.gz | awk '($3=="gene"){print $1":",$4-1000,"-",$5+1000}' | head -500 | tail -1 > ${jbdir}/jbrowse2/plant2tupload/${genome}/defaultRegion.txt
else
    zcat ${wsdir}/${genome}.gff3.gz | awk '($3=="mRNA"){print $1":",$4-1000,"-",$5+1000}' | head -500 | tail -1 > ${jbdir}/jbrowse2/plant2tupload/${genome}/defaultRegion.txt
fi
    
done