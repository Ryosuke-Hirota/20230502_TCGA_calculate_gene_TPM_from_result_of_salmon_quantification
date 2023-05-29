# this script is to convert from transcript TPM to gene TPM 
# 2023/05/02 made

# activate packages for converting
library(tximport)
library(stringr)

# set variables 
# quant.path : directory containing files you want to calculate gene TPM
# Site : site primary which you want to calculate gene TPM
quant.path <-"//fsw-q02/okamura-lab/20221006_TCGA_colon_salmon_quant_transcriptome"
site <-"Colon"

# list transcript quantification files
quant.results <-list.files(path = quant.path,pattern = "quant.sf",recursive = T)

# get id from file path of quantification result
file.name <-str_split(quant.results,pattern = "/",simplify = T)
file.name <-file.name[,1]

# import list of gene id corresponded to transcript id
# this list is located at "https://github.com/Ryosuke-Hirota/20230502_TCGA_calculate_gene_TPM_from_result_of_salmon_quantification"
setwd("C:/Rdata/20230502_TCGA_calculate_gene_TPM_from_result_of_salmon_quantification")
cor.list<-read.table("gencode_v36_list_of_ENST_ID_and_ENSG_ID.txt",sep="\t",header = T,stringsAsFactors = F)

# import list of gene length
# this list is located at "https://github.com/Ryosuke-Hirota/20230502_TCGA_calculate_gene_TPM_from_result_of_salmon_quantification"
length.list<-read.table("gencode_v36_list_of_gene_lenght.txt",sep="\t",header = T,stringsAsFactors = F)

# sum read counts for each gene
setwd(quant.path)
tx.exp <- tximport(quant.results, type = "salmon", txOut = T)
gene.exp <- summarizeToGene(tx.exp,cor.list)
gene.count <-as.data.frame(gene.exp$counts)
gene.count[,ncol(gene.count)+1] <-rownames(gene.count)
rownames(gene.count) <-NULL
gene.count <-gene.count[,c(ncol(gene.count),1:ncol(gene.count)-1)]
colnames(gene.count) <-c("gene_id",file.name)

# calculate TPM for each gene and make table about it 
for (i in 2:ncol(gene.count)) {
  df <-gene.count[,c(1,i)]
  m <-match(df[,1],length.list[,1])
  df[,3] <-length.list[m,2]
  df[,4] <-df[,2]/df[,3]*1000
  df[,5] <-df[,4]/sum(df[,4])*1000000
  colnames(df)[5] <-colnames(df)[2]
  df <-df[,c(1,5)]
  if(i==2){
   tpm.df <-df
    }else{
    tpm.df <-merge(tpm.df,df,by="gene_id")
  }
}
setwd("C:/Rdata/20230502_TCGA_calculate_gene_TPM_from_result_of_salmon_quantification")
write.table(tpm.df,paste0("table_of_gene_TPM_in_",site,".txt"),sep="\t",row.names = F,quote = F)
