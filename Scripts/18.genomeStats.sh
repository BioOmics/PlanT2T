# Usage: bash 18.genomeStats.sh namaList

namaList=$(realpath $1)

db_ID="1"
db_version="v1"
id=$(cat ${namaList} | cut -f1)
name=$(cat metadata.txt | grep "^Latin Name: " | sed 's/Latin Name: //')
common_name=$(cat metadata.txt | grep "^Common Name: " | sed 's/Common Name: //') && if [ -z "${common_name}" ]; then common_name="NA"; fi
full_name=$(cat metadata.txt | grep "^Formatted Species ID: " | sed 's/Formatted Species ID: //')
img_src=$(cat metadata.txt | grep "^PNG File: " | sed 's/PNG File: //g' | sed 's/\/public\/workspace\/biobigdata\/project\/Plant2t/https:\/\/biobigdata.nju.edu.cn\/plant2t/g')
collected_date=$(date '+%Y/%-m/%-d')
Ploidy=$(cat metadata.txt | grep "^Ploidy: " | sed 's/Ploidy: //') && if [ -z "${Ploidy}" ]; then Ploidy="2n"; fi
QV=$(cat metadata.txt | grep "^QV: " | sed 's/QV: //') && if [ -z "${QV}" ]; then QV="NA"; fi
LAI=$(cat metadata.txt | grep "^LAI: " | sed 's/LAI: //') && if [ -z "${LAI}" ]; then LAI="NA"; fi
Email=$(cat metadata.txt | grep "^Email: " | sed 's/Email: //') && if [ -z "${Email}" ]; then Email="NA"; fi
Author=$(cat metadata.txt | grep "^Author: " | sed 's/Author: //') && if [ -z "${Author}" ]; then Author="NA"; fi
Organization=$(cat metadata.txt | grep "^Unit: " | sed 's/Unit: //') && if [ -z "${Organization}" ]; then Organization="NA"; fi
DOI=$(cat metadata.txt | grep "^DOI: " | sed 's/DOI: //') && if [ -z "${DOI}" ]; then DOI="NA"; fi
Message=$(cat metadata.txt | grep "^Message: " | sed 's/Message: //') && if [ -z "${Message}" ]; then Message="NA"; fi
assem_source=$(cat metadata.txt | grep "^Genome Sequencing: " | sed 's/Genome Sequencing: //') && if [ -z "${assem_source}" ]; then assem_source="NA"; fi
GenomeSurvey=$(cat metadata.txt | grep "^Genome Survey: " | sed 's/Genome Survey: //') && if [ -z "${GenomeSurvey}" ]; then GenomeSurvey="NA"; fi
AssemblyByHiFi=$(cat metadata.txt | grep "^Genome Assembly by HiFi: " | sed 's/Genome Assembly by HiFi: //') && if [ -z "${AssemblyByHiFi}" ]; then AssemblyByHiFi="NA"; fi
AssemblyByONT=$(cat metadata.txt | grep "^Genome Assembly by ONT: " | sed 's/Genome Assembly by ONT: //') && if [ -z "${AssemblyByONT}" ]; then AssemblyByONT="NA"; fi
AssemblyChloroplast=$(cat metadata.txt | grep "^Organelle Genomes Assembly: " | sed 's/Organelle Genomes Assembly: //') && if [ -z "${AssemblyChloroplast}" ]; then AssemblyChloroplast="NA"; fi
Polish=$(cat metadata.txt | grep "^Genome Polish: " | sed 's/Genome Polish: //') && if [ -z "${Polish}" ]; then Polish="NA"; fi
HiCScaffolding=$(cat metadata.txt | grep "^Hi-C Scaffolding: " | sed 's/Hi-C Scaffolding: //') && if [ -z "${HiCScaffolding}" ]; then HiCScaffolding="NA"; fi
GapFilling=$(cat metadata.txt | grep "^Gap Filling: " | sed 's/Gap Filling: //') && if [ -z "${GapFilling}" ]; then GapFilling="NA"; fi
QualityAssessment="NA"
TelomereIdentify=$(cat metadata.txt | grep "^Telomere Identification: " | sed 's/Telomere Identification: //') && if [ -z "${TelomereIdentify}" ]; then TelomereIdentify="NA"; fi
CentromereIdentify=$(cat metadata.txt | grep "^Centromere Identification: " | sed 's/Centromere Identification: //') && if [ -z "${CentromereIdentify}" ]; then CentromereIdentify="NA"; fi
Subgenome=$(cat metadata.txt | grep "^Polyploid Subgenome Phasing: " | sed 's/Polyploid Subgenome Phasing: //') && if [ -z "${Subgenome}" ]; then Subgenome="NA"; fi
RepeatAnno=$(cat metadata.txt | grep "^Repeat Annotation: " | sed 's/Repeat Annotation: //') && if [ -z "${RepeatAnno}" ]; then RepeatAnno="NA"; fi
CodingAnno=$(cat metadata.txt | grep "^Gene Model Prediction: " | sed 's/Gene Model Prediction: //') && if [ -z "${CodingAnno}" ]; then CodingAnno="NA"; fi


# 输出标题行
echo -e "db_ID\tdb_version\tid\tname\tcommon_name\tfull_name\timg_src\tcollected_date\tPloidy\ttotal_length\tfragments_num\tN50\tN_count\tGaps\tT2TCount\tT2NCount\tN2NCount\tgene_num\ttranscript_num\tbusco\tTR\tTF\tGC\tEndStop\tMidStop\tNoStop\tncbi_taxonomy_ID\tPhylum\tClass\tOrder\tFamily\tGenus\tQV\tLAI\tEmail\tAuthor\tOrganization\tDOI\tMessage\tassem_source\tGenomeSurvey\tAssemblyByHiFi\tAssemblyByONT\tAssemblyChloroplast\tPolish\tHiCScaffolding\tGapFilling\tQualityAssessment\tTelomereIdentify\tCentromereIdentify\tSubgenome\tRepeatAnno\tCodingAnno"

# 输出数据行
echo -e "${db_ID}\t${db_version}\t${id}\t${name}\t${common_name}\t${full_name}\t${img_src}\t${collected_date}\t${Ploidy}\t$(
    grep -E '^sum|^ave|^largest|^n =|Gaps|^N_count|^N50' genome.stats \
    | sed 's/, /\n/g' \
    | awk -F ' = ' '{print $2}' \
    | tr "\n" "\t" \
    | cut -f1,2,5,7,8
)\t$(
    grep "telomere found" genome.telo.info \
    | cut -d ":" -f2 \
    | cut -d " " -f2 \
    | tr "\n" "\t"
)\t$(
    awk -vFS="\t" '{if($3=="gene"){print $0}}' genome.renamed.gff \
    | wc -l
)\t$(
    awk -vFS="\t" '{if($3=="mRNA"){print $0}}' genome.renamed.gff \
    | wc -l
)\t$(
    cat busco_result/short_summary.*.txt \
    | grep C: \
    | sed -E 's/C:([0-9]*\.[0-9]%).*/\1/g' \
    | cut -f2
)\t$(
    awk -vFS="\t" '{count[$3]++} END {for (type in count) printf "%d\t", count[type]}' genome.re.pep_output/tf_classification.txt
)\t$(
    head -1 genome.re.GC.txt \
    | cut -f2
)\t$(
    cat pep.static.txt
)\t$(
    cat taxonkit.txt \
    | awk -F "[ \t]" -vOFS="\t" '{print $1,$3,$4,$5,$6,$7}'
)\t${QV}\t${LAI}\t${Email}\t${Author}\t${Organization}\t${DOI}\t${Message}\t${assem_source}\t${GenomeSurvey}\t${AssemblyByHiFi}\t${AssemblyByONT}\t${AssemblyChloroplast}\t${Polish}\t${HiCScaffolding}\t${GapFilling}\t${QualityAssessment}\t${TelomereIdentify}\t${CentromereIdentify}\t${Subgenome}\t${RepeatAnno}\t${CodingAnno}" > genome.summary.txt

cat genome.summary.txt
