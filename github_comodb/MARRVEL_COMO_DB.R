library("httr")
library("jsonlite")
library("gprofiler2")

setwd("~/OHSU/conradlab/comodb2")
#mgi<-scan("~/Downloads/mgi_azoospermia_manual_curation_genesonly2.txt",what="character")

#genes<-read.csv("genes/mouse_mgi.txt")


pheno.levels<-c(1:20) #18 MCAT LEVELS PLUS TWO INFERTILITY LEVELS MADE BY US
#APO AND FYPO DON"T HAVE MCAT CATEGORIES?
#FYPO IS POMBE TAXONID 4896
#APO IS CEREVISIAE TAXONID 4932
# 6239 Caenorhabditis elegans
# 7227	Drosophila melanogaster 
# 7955   Danio rerio 
# 10090 Mus musculus
# 10116 Rattus norvegicus
#8364 Xenopus

#FBcv-FlyBase controlled vocabulary (used for phenotype description, allele descrition, etc)
#FBbt-anatomy, not phenotype

#MCAT LEVELS
mcat.lab<-c(1:18)
mcat.lab[1]<-"prenatal_birth"
mcat.lab[2]<-"growth_dev"
mcat.lab[3]<-"nervous"
mcat.lab[4]<-"eye_ear"
mcat.lab[5]<-"integument"
mcat.lab[6]<-"head_neck"
mcat.lab[7]<-"limb"
mcat.lab[8]<-"skeletal_connect"
mcat.lab[9]<-"blood_immune"
mcat.lab[10]<-"cardiovascular"
mcat.lab[11]<-"musculature"
mcat.lab[12]<-"respiratory"
mcat.lab[13]<-"digestive"
mcat.lab[14]<-"genitourinary"
mcat.lab[15]<-"endocrine"
mcat.lab[16]<-"metabolism"
mcat.lab[17]<-"neoplasm"
mcat.lab[18]<-"other"



#LOAD IN SPECIES-SPECIFIC REPRO ANNOTATIONS
flyd<-read.table("fertility_categories/fly_fertility_categories.txt",sep="\t",header=T)
moused<-read.table("fertility_categories/mouse_fertility_categories.txt",sep="\t",header=T)
wormd<-read.table("fertility_categories/worm_fertility_categories.txt",sep="\t",header=T)
zebrad <- read.delim2("~/OHSU/conradlab/comodb2/fertility_categories/zebrafish_fertility_categories.txt")
humand<-read.table("fertility_categories/human_fertility_categories.txt",sep="\t",header=T)

#OBTAIN GENE NAMES TO QUERY
#x<-read.table("genome_mouse.txt",header=T,sep="\t",comment.char="",quote="")
#x<-fread("genome_mouse.txt")
#genes<-unlist(unique(x[,2]))

#genes<-unlist(read.table("genes/mouse_mgi.txt",header=F))
gpro.out<-gconvert(query=as.vector(genes),organism="mmusculus",target="ENTREZGENE_ACC")
entrez.genes<-gpro.out$target

#MAKE OUTPUT CONTAINER
big.out<-NULL

#LOOP OVER GENES
for (k in 1:length(entrez.genes)){
url = paste('http://api.marrvel.org/data/diopt/ortholog/gene/entrezId/',entrez.genes[k],sep="")
res=GET(url)
data = fromJSON(rawToChar(res$content))

if (length(data)==0){next;}

#CONTROL HOMOLOGY RELATIONSHIPS
data<-data[which(data$bestScore==TRUE),]

data<-data[which(data$confidence %in% c("high","moderate")),]

#DROP YEAST FOR NOW, SINCE THERE ARE NO MCAT CATEGORIES
data<-data[! data$taxonId2 %in% c(4896,4932), ]

#PARSE MCAT CATEGORIES AND MAKE OUTPUT 

#FIND COLUMN1
my.col1<-which(names(data)=="gene2")

#FIND COLUMN2
my.column<-which(names(data[[my.col1]])=="phenotypes")
phenos<-data[[my.col1]][my.column]

nspecies<-dim(phenos)[1]

#LOOP OVER EACH SPECIES
#ARE SPECIES ALWAYS SAME ORDER?
big_hit.table<-NULL

#NO PHENOTYPE DATA
if (all(is.na(phenos))){
	big_hit.table<-data.frame(matrix(0,nrow=nspecies,ncol=length(pheno.levels)))
	colnames(big_hit.table)<-as.character(pheno.levels) } else {

#PHENOTYPE DATA
for ( i in 1:nspecies){
  
    
spheno<-phenos[i,][[1]]

if (length(spheno)==0 | length(spheno$ontology)==0){ hit.table<-rep(0,length(pheno.levels));}
else {
#DROP GENE EXPRESSION ANNOTATIONS
    g.hits<-grep("FBbt",spheno$id)
  if(length(g.hits>0)){spheno<-spheno[-g.hits,]}  #DROP FB EXPRESSION ANNOTATION
  spheno<-spheno[,"ontology"]



tmp<-unlist(spheno)
hits<-factor(as.numeric(tmp[tmp %in% pheno.levels]),levels=pheno.levels)
	hit.table<-table(hits)
	
	#SPECIES SPECIFIC ANALYSIS OF FERTILITY PHENOS
	if (data$taxonId2[i]==7227){repro.n<-sum(spheno$id %in% flyd$id); repro.s<-paste(unique(flyd$sum[flyd$id %in% spheno$id]),collapse=";")}
	if (data$taxonId2[i]==10090){repro.n<-sum(spheno$id %in% moused$id); repro.s<-paste(unique(moused$sum[moused$id %in% spheno$id]),collapse=";")}
#	if (data$taxonId2[i]==10116){repro.n<-sum(spheno$id %in% moused$id); repro.s<-paste(unique(moused$sum[moused$id %in% spheno$id]),collapse=";")}
	if (data$taxonId2[i]==6239){repro.n<-sum(spheno$id %in% wormd$id); repro.s<-paste(unique(wormd$sum[wormd$id %in% spheno$id]),collapse=";")}
	if (data$taxonId2[i]==7955){repro.n<-sum(spheno$id %in% zebrad$id); repro.s<-paste(unique(zebrad$sum[zebrad$id %in% spheno$id]),collapse=";")}
	if (data$taxonId2[i]==9606){repro.n<-sum(spheno$id %in% humand$id); repro.s<-paste(unique(humand$sum[humand$id %in% spheno$id]),collapse=";")}
	
	hit.table[19]<-repro.n
	hit.table[20]<-repro.s
}


	big_hit.table<-rbind(big_hit.table,hit.table)
} #END LOOP OVER SPECIES
	  
meta.data<-data[,1:7]
meta.data$gene<-gpro.out$name[k]


	
out<-cbind(meta.data,big_hit.table) }

big.out<-rbind(big.out,out)
}
write.table(big.out,"como_db_v1.txt")


#SIMPLE DOWNSTREAM ANALYSES
#TABULATE RATES ACROSS MCAT CATEGORIES
cats<-apply(big.out[,c(9:26)],MAR=2,FUN=function(x){sum(x>0)})
pdf("MCat counts all species.pdf",width=9,height=6)
barplot(cats,main="All Species All Genes",ylab="Frequency")
dev.off()

x.link<-x$Gene.wgEncodeGencodeBasicV19[which(x$Chr=="X")]
x.cats<-apply(big.out[big.out$gene %in% x.link,c(9:26)],MAR=2,FUN=function(x){sum(x>0)})
barplot(x.cats,main="All Species X-link")
chisq.test(data.frame(c1=cats,c2=x.cats))

#RATES BY SPECIES
m.cats<-apply(big.out[which(big.out$taxonId2 ==10090),c(9:26)],MAR=2,FUN=function(x){sum(x>0)})
r.cats<-apply(big.out[which(big.out$taxonId2 ==10116),c(9:26)],MAR=2,FUN=function(x){sum(x>0)})
z.cats<-apply(big.out[which(big.out$taxonId2 ==7955),c(9:26)],MAR=2,FUN=function(x){sum(x>0)})
w.cats<-apply(big.out[which(big.out$taxonId2 ==6239),c(9:26)],MAR=2,FUN=function(x){sum(x>0)})
f.cats<-apply(big.out[which(big.out$taxonId2 ==7227),c(9:26)],MAR=2,FUN=function(x){sum(x>0)})


#NUMBER OF GENES WITH INFERTILITY PHENOTYPES REPORTED IN EACH SPECIES
ihit<-NULL
ihit[1]<-length(unique((big.out$gene[which(big.out$taxonId2 ==10090 & big.out[,27]>0)])))
ihit[2]<-length(unique((big.out$gene[which(big.out$taxonId2 ==10116 & big.out[,27]>0)])))
ihit[3]<-length(unique((big.out$gene[which(big.out$taxonId2 ==7955 & big.out[,27]>0)])))
ihit[4]<-length(unique((big.out$gene[which(big.out$taxonId2 ==6239 & big.out[,27]>0)])))
ihit[5]<-length(unique((big.out$gene[which(big.out$taxonId2 ==7227 & big.out[,27]>0)])))


########bar plot of gene matches
hits<-NULL
hits[1]<-length(unique(big.out$gene[which(big.out$taxonId2==10090 & big.out$confidence=="high")]))
hits[2]<-length(unique(big.out$gene[which(big.out$taxonId2==10116 & big.out$confidence=="high")]))
hits[3]<-length(unique(big.out$gene[which(big.out$taxonId2==7955 & big.out$confidence=="high" )]))
hits[4]<-length(unique(big.out$gene[which(big.out$taxonId2==6239 & big.out$confidence=="high")]))
hits[5]<-length(unique(big.out$gene[which(big.out$taxonId2==7227& big.out$confidence=="high" )]))
hits[6]<-length(unique(big.out$gene[which(big.out$taxonId2==8364 & big.out$confidence=="high")]))

barplot(hits)