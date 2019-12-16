
# Author: Annabel Smith

plink_locus_info<-function(locus_info_data){

# locus_data = full locus info data frame from STEP_01

# LOCUS INFO:
linf2<-data.frame(locus_info_data[,c("locus","AlleleID")])
head(linf2)

# Add allele columns:
linf2$ref_allele<-substr(linf2$AlleleID,nchar(as.character(linf2$AlleleID))-2,nchar(as.character(linf2$AlleleID))-2)
linf2$snp_allele<-substr(linf2$AlleleID,nchar(as.character(linf2$AlleleID)),nchar(as.character(linf2$AlleleID)))
check.rows(linf2)

# Add genotypes to locus info:
# 0 = Reference allele homozygote
# 1 = SNP allele homozygote
# 2 = heterozygote
linf2$gt0<-paste(linf2$ref_allele,linf2$ref_allele)
linf2$gt1<-paste(linf2$snp_allele,linf2$snp_allele)
linf2$gt2<-paste(linf2$ref_allele,linf2$snp_allele)

head(linf2)
return(linf2)

} # close plink locus info

format_plink_ped<-function(snp_data,locus_data,remove_og=NULL,remove_cultivar=NULL){

# snp_data = imported and formatted snp data (from STEP_01)
# locus_data = linf, all locus info from STEP_01
# remove_og = remove outgroups?
# remove_cultivar = remove culivars?

print(paste("Start time = ",Sys.time(),sep=""))

# ALLELE ORDER FOR HETEROZYGOTES:
# "allele order in a .ped heterozygous genotype does not matter, and PLINK does not keep track of it" https://groups.google.com/forum/#!searchin/plink2-users/order$20heterozygous$20allele|sort:relevance/plink2-users/eb-XObwwkJU/veWukEB_9wsJ
# "plink never tries to maintain heterozygous call allele order in .ped files; "A B" is considered equivalent to "B A" and they're interchanged arbitrarily.  It's necessary to use a different data representation if order is important" https://groups.google.com/forum/#!searchin/plink2-users/order$20heterozygous$20allele|sort:relevance/plink2-users/IWRMZY8XeyM/k8b2NlIrDAAJ

## ~~~~ ****** .ped file ****** ~~~~ ##

# The .ped file can be can be SPACE or TAB delimited: http://www.shapeit.fr/pages/m02_formats/pedmap.html
# https://www.cog-genomics.org/plink2/formats
# https://www.researchgate.net/post/How_do_I_convert_a_SNP_genotype_table_into_plink_binary_PED_files

# run plink locus info function
linf2<-plink_locus_info(locus_data)

# Use single row data, mono removed:
ghead(snp_data)

# data.frame for the output:
plink1<-snp_data
all_loci<-colnames(plink1)[3:length(plink1)]
plink.list1<-list()

# Convert genotypes:
# ~ 11 MINS for full data set (DPlan18)

for (i in 1:length(all_loci)){

locus.thisrun<-all_loci[i]

data.thisrun<-plink1[,locus.thisrun]
linf.thisrun<-linf2[linf2$locus==locus.thisrun,]

plink.list1[[i]]<-ifelse(is.na(ifelse(ifelse(ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun)==1,linf.thisrun$gt1,ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun))==2,linf.thisrun$gt2,ifelse(ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun)==1,linf.thisrun$gt1,ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun)))),"0 0",ifelse(ifelse(ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun)==1,linf.thisrun$gt1,ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun))==2,linf.thisrun$gt2,ifelse(ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun)==1,linf.thisrun$gt1,ifelse(data.thisrun==0,linf.thisrun$gt0,data.thisrun))))

} # close for loci

plink2<-data.frame(plink1[,1:2],do.call(cbind,plink.list1))
colnames(plink2)[3:length(plink2)]<-all_loci
ghead(plink2)

# Add parents, sex & phenotype (0 = unknown):
plink2$PaternalID<-rep(0,nrow(plink2))
plink2$MaternalID<-rep(0,nrow(plink2))
plink2$Sex<-rep(0,nrow(plink2))
plink2$Phenotype<-rep(0,nrow(plink2))
ghead(plink2)

# Re-arrange:
plink2<-plink2[,c(which(colnames(plink2) %in% c("site","ind","PaternalID","MaternalID","Sex","Phenotype")),head(grep("L",colnames(plink2)),1):tail(grep("L",colnames(plink2)),1))]

if(is.null(remove_og)==F){
if (remove_og==T){
if(length(grep("OG",plink2$site))>0){
# remove outgroups:
plink2<-plink2[-grep("OG",plink2$site),]
plink2<-tidy.df(plink2)
} # close if og
} # close remove_og
} # close is null

if(is.null(remove_cultivar)==F){
if (remove_cultivar==T){
cultivar_rows<-c(grep("CAT",plink2$site),grep("CCT",plink2$site),grep("CTP",plink2$site))
if(length(cultivar_rows)>0){
# remove culitvar:
plink2<-plink2[-cultivar_rows,]
plink2<-tidy.df(plink2)
} # close if culivar
} # close remove_cultivar
} # close is null

print(paste("End time = ",Sys.time(),sep=""))

return(plink2)

} # close plink .ped

check_plink_ped<-function(orig_data,plink_data){

# orig_data = the inputted, formatted snp data from STEP_01
# plink_data = plink formatted ped file from format_plink_ped function

# check:
locus.now<-sample(colnames(orig_data)[3:length(orig_data)],3)
ind.now<-sample(unique(orig_data$ind),3)
print(orig_data[ind.now,c(1,2,which(colnames(orig_data) %in% locus.now))][order(orig_data[ind.now,c(1,2,which(colnames(orig_data) %in% locus.now))]$ind),])
print(plink_data[which(plink_data$ind %in% ind.now),c(1,2,which(colnames(plink_data) %in% locus.now))][order(plink_data[which(plink_data$ind %in% ind.now),c(1,2,which(colnames(plink_data) %in% locus.now))]$ind),])
# 0 = Reference allele homozygote
# 1 = SNP allele homozygote
# 2 = heterozygote

} # close check plink ped

format_plink_map<-function(ped_file,locus_data){

# ped_file = ped file formatted with the format_plink_ped function
# locus_data = linf, all locus info from STEP_01

# Subset locus data using loci in the ped_file:
loci_now<-colnames(ped_file)[grep("L",colnames(ped_file))[1]:length(ped_file)]
linf_now<-locus_data[locus_data$locus %in% loci_now,]
linf_now<-tidy.df(linf_now)

# Should be zero
if(length(which(!colnames(ped_file)[grep("L",colnames(ped_file))] %in% linf_now$locus))!=0) stop("Something is wrong")

map_now<-data.frame(chrom=rep(1,nrow(linf_now)),snp=linf_now$locus,distance=rep(0,nrow(linf_now)),position=linf_now$SnpPosition)

return(map_now)

} # close plink map

write_parameters<-function(data_name,sites_now,description,out_file)
{
cat(description,
paste("Time = ",Sys.time(),sep=""),
paste("Data name = ",data_name,sep=""),
paste("Number of sites = ",length(sites_now),sep=""),
paste("Sites included = ",paste(sites_now,collapse=", "),sep=""),
paste("Number of individuals = ",nrow(get(data_name)),sep=""),
paste("Monomorphic loci removed = ",remove_mono,sep=""),
paste("Maximum missing data = ",cr*100,"%",sep=""),
paste("Minimum reproducibility = ",ra*100,"%",sep=""),
paste("Minimum minor allele frequency = ",malim*100,"%",sep=""),
file=paste(out_file,".txt",sep=""),sep="\n",append=T)

} # close write_parameters

read_prune_out<-function(directory,correlation){
LD_df<-as.character(unlist(read.table(paste(directory,dir(directory)[grep(paste(correlation,".prune.out",sep=""),dir(directory))],sep="/"),header=F)))
return(LD_df)
}

