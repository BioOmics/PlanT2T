# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否传递了参数
if (length(args) == 0) {
  stop("Usage: Rscript OrgDBmaker.R Arabidopsis_thaliana_Col-CEN")
}

# install.packages("")
library(dplyr)
library(stringr)
library(jsonlite)
library(AnnotationForge)
options(stringsAsFactors = F)


# 传递的物种名称
speciesID <- args[1]
species <- args[2]
tax_id <- args[3]


# 构建文件名
file_name <- paste0(speciesID, "_gene.csv")
emapper <- read.table(file_name, header=TRUE, sep = ",", quote = "")
emapper[emapper==""]<-NA

gene_info <- emapper %>% dplyr::select(GID = gene_ID, GENENAME = source_ID) %>% na.omit()
gos <- emapper %>% dplyr::select(gene_ID, go) %>% na.omit()

gene2go = data.frame(GID = character(),
                     GO = character(),
                     EVIDENCE = character())
                     
gos_list <- function(x){
  the_gos <- str_split(x[2], ";", simplify = FALSE)[[1]]
  df_temp <- data.frame(GID = rep(x[1], length(the_gos)),
                        GO = the_gos,
                        EVIDENCE = rep("IEA", length(the_gos)))
  return(df_temp)
}

gene2gol <- apply(as.matrix(gos),1,gos_list)
gene2gol_df <- do.call(rbind.data.frame, gene2gol)
gene2go <- gene2gol_df
gene2go$GO[gene2go$GO=="-"]<-NA
gene2go<-na.omit(gene2go)


gene2ko <- emapper %>% dplyr::select(GID = gene_ID, Ko = kegg)
gene2ko$Ko[gene2ko$Ko=="-"]<-NA
gene2ko<-na.omit(gene2ko)
gene2kol <- apply(as.matrix(gene2ko),1,gos_list)
gene2kol_df <- do.call(rbind.data.frame, gene2kol)
gene2ko <- gene2kol_df[,1:2]
colnames(gene2ko) <- c("GID","Ko")
gene2ko$Ko <- gsub("ko:","",gene2ko$Ko)

ko2pathway <- read.table("11.ko00001.filter4Plant.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", check.names = F)[, c(5, 4)]
colnames(ko2pathway) <- c("Ko", "Pathway")

gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% dplyr::select(GID, Pathway) %>% na.omit()

genus = "" 
species = gsub("_", "\\.", species)
species = gsub("-", "\\.", species)

gene2go <- unique(gene2go)
gene2go <- gene2go[!duplicated(gene2go),]
gene2ko <- gene2ko[!duplicated(gene2ko),]
gene2pathway <- gene2pathway[!duplicated(gene2pathway),]
gene_info <- gene_info[!duplicated(gene_info),]

makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               pathway=gene2pathway,
               version="1.0",
               maintainer = "PlanT2T <haoyuchao@zju.edu.cn>",
               author = "PlanT2T",
               outputDir = ".",
               tax_id=tax_id,
               genus=genus,
               species=species,
               goTable="go")
