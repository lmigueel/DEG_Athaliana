# Desenvolvedor: Lucas Miguel  de Carvalho - UNICAMP #
# lucasmiguel@lge.ibi.unicamp.br #
# Script de analise de dados do artigo: #
# https://www.ncbi.nlm.nih.gov/bioproject/555093 
# Acessar os dados do SRA ###
# https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=555093 #


# Primeiro passo seria baixar o arquivo SraRunInfo.tsv#
# Ele contem todos os links do SRA eo ID das amostras #
# Selecionamos as mostras SRR9696658, SRR9696662, SRR9696666,SRR9696660,SRR9696664,SRR9696668
# posteriormente clicar em 'Send to' -> File -> RunInfo 

###### SRA #######

setwd(".")
#https://www.ncbi.nlm.nih.gov/sra/?term=SRP215218
base_dir <- getwd()

dados <-read.csv("SraRunInfo.csv", stringsAsFactors=FALSE)
arquivos <- basename(dados$download_path)
for(i in 1:length(arquivos)){
  download.file(dados$download_path[i], arquivos[i])
}


for(a in arquivos) {
  cmd = paste("fastq-dump --split-3", a)
  system(cmd)
}


######  Trimmomatic #####
#http://www.usadellab.org/cms/?page=trimmomatic

cmd = paste("wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip")
system(cmd)

cmd = paste("unzip Trimmomatic-0.39.zip")
system(cmd)

for(a in arquivos){ 
  cmd = paste("java -jar ",base_dir,"/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -trimlog ",a,".trimlog -summary ",a,".summary ",a,".fastq ",a,".trim.fastq ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",sep = "")
  system(cmd)
}

######  FastQC ######
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

cmd = paste("wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip")

cmd = paste("unzip fastqc_v0.11.8.zip")
system(cmd)

cmd = paste("chmod 755 ",base_dir,"FastQC/fastqc")
system(cmd)

for(a in arquivos){ 
  cmd = paste(base_dir,"FastQC/fastqc ",a,".fastq ",sep = "")
  system(cmd)
}

###### Kallisto ####

cmd = paste("wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz")

cmd = paste("tar -xzvf kallisto_linux-v0.44.0.tar.gz")
system(cmd)

## Voce deve criar um indice com o kallisto de seu transcriptoma
## comando:  kallisto index -i arabidopsis_index <transcriptoma>
for(a in arquivos){ 
  cmd = paste(base_dir,"/kallisto_linux-v0.44.0/kallisto quant -i arabidopsis_index -o ",a,"_kallisto -b 100 -t 10 --single -l 100 -s 0.001 ",a,".trim.fastq",sep="")
  system(cmd)
}

##### Copiar arquivos ########

cmd = paste("for file in ls -l -d SRR*;do cp $file/abundance.tsv $file.tsv;done")
print(cmd)
system(cmd)

###### Montar matriz ########

# precisa de um script presente no trinity 
# cmd = ("wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.8.6/trinityrnaseq-v2.8.6.FULL.tar.gz")
# system (cmd)

lista = paste0(arquivos,".tsv",collapse = " ")
cmd = paste("perl trinityrnaseq-2.8.6/util/abundance_estimates_to_matrix.pl  --est_method kallisto --gene_trans_map none ",lista,sep="")
print (cmd)

#gera uma matriz chamada kallisto.isoform.counts.matrix que sera utilizada
#no deseq2 e edgeR

####### Sleuth #########

#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#install.packages("devtools", repos = "http://cran.us.r-project.org")
#library("httr")
#set_config(config(ssl_verifypeer = 0L))
#devtools::install_github("pachterlab/sleuth")
library("sleuth")

cmd = "mkdir sleuth"
system(cmd)
cmd = "cp -r SRR*/ sleuth/"
system(cmd)

base_dir <- getwd()

########## SLEUTH #############
#HS1     SRR9696660
#HS2     SRR9696664
#HS3     SRR9696668
#CT1     SRR9696658      
#CT2     SRR9696662
#CT3     SRR9696666

sample_id <- list('SRR9696658','SRR9696662','SRR9696666',
                  'SRR9696660','SRR9696664','SRR9696668')

paths <- list(paste(base_dir,"/sleuth/SRR9696658",sep=""),
              paste(base_dir,"/sleuth/SRR9696662",sep=""),
              paste(base_dir,"/sleuth/SRR9696666",sep=""),
              paste(base_dir,"/sleuth/SRR9696660",sep=""),
              paste(base_dir,"/sleuth/SRR9696664",sep=""),
              paste(base_dir,"/sleuth/SRR9696668",sep=""))

names(paths) <- sample_id

s2c <- read.table(file.path(base_dir, "amostras.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = sample, condition, reps)
s2c
#t2g <- read.table("t2g.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = paths)
print(s2c)

s2c <- data.frame(lapply(s2c, as.character), stringsAsFactors=FALSE)

#transcrito
so <- sleuth_prep(s2c, ~condition, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so)

#wald
so <- sleuth_wt(so, "conditionTratado")
models(so)
results_table <- sleuth_results(so, test='conditionTratado', test_type = 'wald')
sleuth_significant <- dplyr::filter(results_table, qval <= 0.05)
head(sleuth_significant, 20)
sleuth_list <- sleuth_significant[,1]


write.table(sleuth_significant,file="diferenciais_sleuth.txt")

sleuth_live(so)

pdf("Sleuth_Volcano.pdf")
plot(results_table$b, -1*log10(results_table$qval), col=ifelse(results_table$qval, "red", "black"),xlab="log(qval)", ylab="beta", title="Volcano plot", pch=20)
dev.off()

##### EDEGR #######

if (! require(edgeR)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  library(edgeR)
}

data = read.table("kallisto.isoform.counts.matrix", header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = factor(c(rep("Controle", 3), rep("Tratado", 3)))

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateDisp(exp_study)
et = exactTest(exp_study, pair=c("Controle", "Tratado"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="Controle", sampleB="Tratado", result_table)
result_table$logFC = -1 * result_table$logFC
write.table(result_table, file='edgeR.DE_results', sep='     ', quote=F, row.names=T)
write.table(rnaseqMatrix, file='edgeR.count_matrix', sep='   ', quote=F, row.names=T)

pdf("edgeR.Volcano.pdf")
plot(result_table$logFC, -1*log10(result_table$FDR), col=ifelse(result_table$FDR<=0.05, "red", "black"),xlab="logCounts", ylab="logFC", title="Volcano plot", pch=20)
dev.off()

edger_significant <- dplyr::filter(result_table, FDR <= 0.05)
edgeR_list <- row.names(edger_significant)[i]

edgeR_list <- NULL
for(i in 1:length(result_table[,1])){
    if(result_table[i,]$FDR <= 0.05){
      #print("Entrou\n")
      edgeR_list[i] <- row.names(result_table)[i]
  }
}

head(edgeR_list)

############### DESEQ2 ##########

if (! require(DESeq2)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
}

data = read.table("kallisto.isoform.counts.matrix", header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("Controle", 3), rep("Tratado", 3))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","Controle","Tratado")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Controle"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Tratado"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Controle", sampleB="Tratado", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='DESeq2.DE_results', sep='     ', quote=FALSE)

pdf("DESeq2_Volcano.pdf")
plot(res$log2FoldChange, -1*log10(res$padj), col=ifelse(res$padj<=0.05, "red", "black"),xlab="logCounts", ylab="logFC", title="Volcano plot", pch=20)
dev.off()

Deseq2_list <- NULL
for(i in 1:length(res[,1])){
  if(res[i,]$padj <= 0.05){
    #print("Entrou\n")
    Deseq2_list[i] <- row.names(res)[i]
  }
}

head(Deseq2_list)

#######  Volcano Plot #########


gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE)
layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 
plot(res$log2FoldChange, -1*log10(res$padj), col=ifelse(res$padj<=0.05, "red", "black"),xlab="logCounts", ylab="logFC", title="Volcano plot", pch=20)
plot(result_table$logFC, -1*log10(result_table$FDR), col=ifelse(result_table$FDR<=0.05, "red", "black"),xlab="logCounts", ylab="logFC", title="Volcano plot", pch=20)
plot(sleuth_significant$b, -1*log10(sleuth_significant$qval), col=ifelse(sleuth_significant$qval, "red", "black"),xlab="log(qval)", ylab="beta", title="Volcano plot", pch=20)

am <- read.table(file="amostras.txt",sep="\t",header=TRUE)
head(am)

########### Diagrama de Venn #######

sleuth_significant <- read.table(file="diferenciais_comp1.txt",sep=" ",header=T)
sleuth_list <- sleuth_significant[,1]

edgeR_list <- read.table(file="lista_diff_edgeR.txt")

install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Paired")


venn.diagram(
  x = list(edgeR_list,Deseq2_list,sleuth_list),
  category.names = c("edgeR" , "Deseq2","Sleuth"),
  filename = 'venn_diagramm_DEG.png',
  output=FALSE,
  
  #Saida
  imagetype="png" ,
  height = 680 , 
  width = 880 , 
  resolution = 600,
  
  #Numeros
  cex = .5,
  fontface = "bold",
  fontfamily = "serif",

  #Circulos
  lwd = 2,
  lty = 'blank',
  fill = myCol,

  
  #Nomes
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "serif",
  rotation = 1
  
)

#########  PCA #######

library("DESeq")
countsTable <- read.delim("kallisto.isoform.counts.matrix", header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable[,1]
countsTable <- countsTable[,2:7]
conds <- factor(c(names(countsTable)))

countsTable_novo <- apply(countsTable,2,as.integer)
countsTable_novo[is.na(countsTable_novo)] <- 0

cds<-newCountDataSet(countsTable_novo,conds)
cds<-estimateSizeFactors(cds)
sizeFactors(cds)
cds <- estimateDispersions(cds,method='blind')
vsd <- varianceStabilizingTransformation(cds)

pdf("PCA.pdf")
plotPCA(vsd)
dev.off()

######### DENDOGRAMA #########

install.packages("ggdendro")
install.packages('dendextend') 
library('dendextend')
library("ggplot2")
library("ggdendro")

countsTable <- read.delim("kallisto.isoform.counts.matrix", header=TRUE, stringsAsFactors=TRUE,row.names = 1)
dd <- dist(t(scale(countsTable)), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE,
             size  = 1) + labs(title="Dendrogram in ggplot2")+
              xlab("Amostras") +ylab("Altura")

