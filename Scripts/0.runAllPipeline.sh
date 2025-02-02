#!/bin/bash
ulimit -s 10240000

GenomePathway=$(realpath "$1")

function CheckSoftware() {
    if command -v "$1" >/dev/null 2>&1; then
        sleep 0.05
    else
        echo -e "Error: ${2} can not be found in the PATH environment variable."
        exit 1
    fi
}

CheckSoftware "faSize" "faSize"
CheckSoftware "gffread" "gffread"
CheckSoftware "STAR" "STAR"
CheckSoftware "rsem-prepare-reference" "rsem"
CheckSoftware "faToTwoBit" "faToTwoBit"
CheckSoftware "assembly-stats" "assembly-stats"
CheckSoftware "exec_annotation" "kofamscan"
CheckSoftware "busco" "busco"
CheckSoftware "Rscript" "Rscript"
CheckSoftware "R" "R"
CheckSoftware "python" "python"
CheckSoftware "perl" "perl"
CheckSoftware "quartet.py" "quarTeT"
CheckSoftware "taxonkit" "taxonkit"
CheckSoftware "makeblastdb" "makeblastdb"
CheckSoftware "samtools" "samtools"
CheckSoftware "mysql" "mysql"


cd "$GenomePathway"

# Show metadata.txt
echo "Your form"
echo "----------------------------------"
cat metadata.txt | sed 's/\/public\/workspace\/biobigdata\/project\/Plant2t\/UserUpload//g'

awk -F': ' '
BEGIN { OFS="\t" }
$1 ~ /^ID/ { id=$2 }
$1 ~ /^Formatted Name/ { formatted_name=$2 }
$1 ~ /^Formatted Species ID/ { FormattedSpeciesID=$2 }
END {
    print id, formatted_name, FormattedSpeciesID
}
' metadata.txt > namelist.txt 

namaList=namelist.txt

if [ ! -f "${namaList}" ] || [ ! -s "${namaList}" ]; then
    echo "File not found or empty: ${namaList}"
    exit 1
fi

echo "------------------------------------------"
echo "ID    GeneNamePrefix FormattedSpeciesID" && cat namelist.txt

# Auto increment step number
step=0
function StepCounter() {
    echo -e "------------------------------------------"
    echo -e "Step${step}: $1"
    echo -e "------------------------------------------"
    ((step++))
}

# Wrap each step in a function call
function RunStep() {
    StepCounter "$1"
    bash "$2" "$namaList"
    if [ $? -ne 0 ]; then
        echo -e "Error: $1 failed."
        # rm -rf $GenomePathway
        exit 1
    fi
}

RunStep "00.genomeCheck" "00.genomeCheck.sh"
RunStep "01.makePEPCDS" "01.makePEPCDS.sh"
RunStep "02.renameGff" "02.renameGff.sh"
RunStep "03.teloExplorer" "03.teloExplorer.sh"
RunStep "04.rsemIndex" "04.rsemIndex.sh"
RunStep "05.genome2bit" "05.genome2bit.sh"
RunStep "06.assemblyStats" "06.assemblyStats.sh"
RunStep "07.tfIdent" "07.tfIdent.sh"
RunStep "08.runInterProScan" "08.runInterProScan.sh"
RunStep "09.runKoFamScan" "09.runKoFamScan.sh"
RunStep "10.runBUSCO" "10.runBUSCO.sh"
RunStep "11.orgDBmaker" "11.orgDBmaker.sh"
RunStep "12.centroMiner" "12.centroMiner.sh"
RunStep "13.ideogram" "13.ideogram.sh"
RunStep "14.getKEGG" "14.getKEGG.sh"
RunStep "15.cleanPEP" "15.cleanPEP.sh"
RunStep "16.pepStatic" "16.pepStatic.sh"
RunStep "17.taxonkitFinder" "17.taxonkitFinder.sh"
RunStep "18.genomeStats" "18.genomeStats.sh"
RunStep "19.protparamAnalysis" "19.protparamAnalysis.sh"
RunStep "20.SaveResultToMySQL" "20.SaveResultToMySQL.sh"
RunStep "21.makeBlastDB" "21.makeBlastDB.sh"
RunStep "22.JBrowse2" "22.JBrowse2.sh"
RunStep "23.downloadFile" "23.downloadFile.sh"

echo -e "All done!"
