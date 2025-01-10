#!/bin/bash

read -r -d '' CODE_TO_ADD << EOM
    ## Load pathway data 'ko00001.PlanT2T.txt' as 'pathway4plant'
    data_file <- system.file("extdata", "ko00001.PlanT2T.txt", package = pkgname)
    if (file.exists(data_file)) {
        pathway4plant <- read.table(data_file, header = TRUE, sep = "takeplacename", stringsAsFactors = F, quote = "", check.names = F)
        assign("pathway4plant", pathway4plant, envir = .GlobalEnv)
        packageStartupMessage("Loaded 'pathway4plant' from ko00001.PlanT2T.txt for KEGG annotation")
		packageStartupMessage(" ")
		packageStartupMessage("See more information in PlanT2T:")
		packageStartupMessage("NJU site: https://biobigdata.nju.edu.cn/plant2t/")
		packageStartupMessage("ZJU site: https://bis.zju.edu.cn/plant2t/")
    } else {
        packageStartupMessage("Pathway data 'ko00001.PlanT2T.txt' not found.")
		packageStartupMessage("If you have any questions, please open an issue on GitHub")
		packageStartupMessage("GitHub site: https://github.com/BioOmics/PlanT2T")
    }
EOM

OrgDB=${1}
ZZZ_R_FILE="${OrgDB}/R/zzz.R"

if [[ ! -f "$ZZZ_R_FILE" ]]; then
    echo "Error: $ZZZ_R_FILE does not exist."
    exit 1
fi

onLoadStartLine=$(grep -n '.onLoad <- function' "$ZZZ_R_FILE" | cut -d ':' -f 1)
onLoadBraceLine=$(awk "/^ *{/{print NR; exit}" "$ZZZ_R_FILE")

if [[ -z "$onLoadStartLine" ]] || [[ -z "$onLoadBraceLine" ]]; then
    echo "Error: .onLoad function or its opening brace not found in $ZZZ_R_FILE."
    exit 1
fi

if grep -q "Load pathway data 'ko00001.PlanT2T.txt'" "$ZZZ_R_FILE"; then
    echo "The code block already exists in $ZZZ_R_FILE. Skipping insertion."
else
    awk -v code="$CODE_TO_ADD" -v braceLine="$onLoadBraceLine" '
    NR == braceLine { print; print code; next }
    { print }
    ' "$ZZZ_R_FILE" > "${ZZZ_R_FILE}.tmp" && mv "${ZZZ_R_FILE}.tmp" "$ZZZ_R_FILE"
    
    echo "Code block inserted into $ZZZ_R_FILE."

    sed -i 's/takeplacename/\\t/g' "$ZZZ_R_FILE"
fi

