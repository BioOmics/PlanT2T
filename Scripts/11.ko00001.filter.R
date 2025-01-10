# Usage: /bioapps/rhel7/R-4.0.0/bin/Rscript 11.ko00001.filter.R -i 11.ko00001.json
package_list <- c("optparse","rjson","stringr","tidyr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#clean enviroment object
rm(list=ls()) 

#Load essential packages
library(optparse)
library(rjson)
library(stringr)
library(tidyr)

option_list = list(make_option(c("-i", "--input"), type = "character", default = NULL,
                               help="Rscript 11.ko00001.filter.R -i 11.ko00001.json"))         
apt = parse_args(OptionParser(option_list=option_list))

jsonData = fromJSON(file = apt$input)
a=jsonData$name
b=jsonData$children

c=NULL
for (i in 1:length(b)) {
  for (j in 1:length(b[[i]]$children)) {
    for (w in 1:length(b[[i]]$children[[j]]$children)) {
      for (x in 1:length(b[[i]]$children[[j]]$children[[w]]$children)) {
        
        c=rbind(c,c(
          b[[i]]$name,
          b[[i]]$children[[j]]$name,
          b[[i]]$children[[j]]$children[[w]]$name,
          b[[i]]$children[[j]]$children[[w]]$children[[x]]$name
        ))
        
      }
    }  
  }
}

c=as.data.frame(c)
c$V1=substring(c$V1,7,nchar(c$V1))
c$V2=substring(c$V2,7,nchar(c$V2))
c$V3=substring(c$V3,7,nchar(c$V3))

d = separate(data=c, col = V3,  into=c('V3', 'pathway_id'), sep ="[[]" )
d$pathway_id <- str_replace_all(d$pathway_id, c("PATH:" = "", "BR:" = "", "]" = ""))

d$V4 =sub("  ", "$", d$V4)
d = separate(data=d, col = V4,  into=c('V4', 'KO_name'), sep ="[$]" )

d$KO_name =sub(";", "$", d$KO_name)
d = separate(data=d, col = KO_name,  into=c('gene', 'KO_name'), sep ="[$]" )

d$KO_name =sub("[[]EC:", "$", d$KO_name)
d = separate(data=d, col = KO_name,  into=c('KO_name', 'EC'), sep ="[$]" )
d$EC =gsub("[]]", "", d$EC)

names(d)[c(1,2,3,5)]=c("level1","level2","level3","KO")

d$KO_name <- trimws(d$KO_name)
d$level3 <- sub("\\s+$", "", d$level3)

d <- subset(d, !(level1 %in% c("Human Diseases", "Organismal Systems")))
d <- subset(d, !(level2 %in% c("Signaling molecules and interaction","Cellular community - eukaryotes","Information processing in viruses", "Cellular community - prokaryotes", "Not included in regular maps", "Unclassified: metabolism", "Unclassified: genetic information processing", "Unclassified: signaling and cellular processes", "Poorly characterized")))
exclude_paths <- c(
  "ErbB signaling pathway", 
  "Ras signaling pathway", 
  "Rap1 signaling pathway", 
  "Wnt signaling pathway", 
  "Notch signaling pathway", 
  "Hedgehog signaling pathway", 
  "TGF-beta signaling pathway", 
  "Hippo signaling pathway", 
  "VEGF signaling pathway", 
  "Apelin signaling pathway", 
  "JAK-STAT signaling pathway", 
  "NF-kappa B signaling pathway", 
  "TNF signaling pathway", 
  "HIF-1 signaling pathway", 
  "FoxO signaling pathway",
  "p53 signaling pathway",
  "Cytoskeleton in muscle cells",
  "Regulation of actin cytoskeleton",
  "Flagellar assembly",
  "Cilium and associated proteins",
  "Domain-containing proteins not elsewhere classified",
  "CD molecules",
  "Prokaryotic defense system",
  "Fanconi anemia pathway",
  "Cytokines and neuropeptides",
  "Necroptosis",
  "Oocyte meiosis",
  "Efferocytosis",
  "Phosphotransferase system (PTS)",
  "cAMP signaling pathway",
  "Phagosome",
  "Lipopolysaccharide biosynthesis proteins",
  "Peptidoglycan biosynthesis and degradation proteins",
  "Antimicrobial resistance genes",
  "Steroid hormone biosynthesis",
  "Arachidonic acid metabolism",
  "Mycolic acid biosynthesis"
)
d <- subset(d, !(level3 %in% exclude_paths))
d <- d[!grepl("- animal|- yeast|- other|- Caulobacter|- fly|Bacterial|- multiple species|Insect|Viral", d$level3), ]

d <- d[!is.na(d$KO), ]

ko2K2name <- d[, c("pathway_id", "KO", "level3")]
colnames(ko2K2name) <- c("Pathway", "Ko", "Name")

write.table(d, file = "11.ko00001.filter4Plant.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(ko2K2name, "14.ko2K2name.txt", quote = F, sep = "\t", row.names = F, col.names = T)
