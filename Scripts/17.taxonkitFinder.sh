# Usage: bash 17.taxonkitFinder.sh namaList
i=$(cat metadata.txt | grep "^Latin Name: " | sed 's/Latin Name: //')

echo ${i} | awk -F "[ _\t]" -vOFS=" " '{print $1,$2}' \
| taxonkit name2taxid | awk -F"\t" '$2!=""{print $2}' \
| taxonkit lineage | taxonkit reformat -F | cut -f 1,3 | sed 's/;/\t/g' > taxonkit.txt

# 如果为空，补NA
if [ ! -s "taxonkit.txt" ]; then
    echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA" > taxonkit.txt
else
    cat taxonkit.txt
fi

cat taxonkit.txt
