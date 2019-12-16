
# Author: Annabel Smith, unless other source indicated

# filter monomorphic loci:
mono_loci<-function(data,first_col){

# data = a dataframe of dart SNPs, cleaned and formatted in STEP_01
# first_col = the first column containing SNP data, i.e. not site or individual data

gt_summary<-data.frame(locus=colnames(data[,first_col:length(data)]),no.genotypes=apply(data[,first_col:length(data)],2,function(x)length(unique(x[!is.na(x)]))))
mono_f1<-as.character(gt_summary$locus[gt_summary$no.genotypes==1])
filt1_mono<-data[,-which(colnames(data) %in% mono_f1)]
filt1_mono<-tidy.df(filt1_mono)

print(paste("Original data: ",nrow(data)," individuals; ",length(data)-2," loci",sep=""))
print(paste("Filtered data: ",nrow(filt1_mono)," individuals; ",length(filt1_mono)-2," loci",sep=""))
print(paste((length(data)-2)-(length(filt1_mono)-2)," loci monomorphic removed",sep=""))

return(filt1_mono)

} # close mono_loci

missing_data<-function(data,first_col,max_missing){

# data = a dataframe of dart SNPs, cleaned and formatted in STEP_01
# first_col = the first column containing SNP data, i.e. not site or individual data
# min_allowed = the minimum allowed missing data proportion (between 0:1) 

# Remarks:
# Cannot use DartSeq CallRate filter because it was calculated on the whole dataset. Individuals from the RCH and FS populations will cause higher proportion of missing data in the rest of the data:
# This metric is highly correlated with the DartSeq metric (99%), but it filters loci biased by the presence of other populations, e.g. removing the firescape samples results in approx. 500 extra loci filtered.
# cr<-linf[,c("locus","CallRate")]
# miss<-merge(miss,cr,by="locus",all.x=T,all.y=F)
# plot(miss$missing_data,miss$CallRate)
# cor.test(miss$missing_data,miss$CallRate)
# callrate50<-miss$locus[which(miss$missing_data>0.5)]
# callrate50_dart<-miss$locus[which(linf$CallRate<0.5)]

md_df<-apply(data[,first_col:length(data)],2,function(x) length(which(is.na(x)))/length(x))
md_df<-data.frame(locus=names(md_df),missing_data=md_df)
md_df<-tidy.df(md_df)

loci_to_filter<-md_df$locus[which(md_df$missing_data>max_missing)]

if(length(loci_to_filter)==0){
print(paste("Original data: ",nrow(data)," individuals; ",length(data)-2," loci",sep=""))
print(paste("Filtered data: ",nrow(data)," individuals; ",length(data)-2," loci",sep=""))
print(paste("There were no loci with more than ", max_missing*100," % missing data; none removed",sep=""))
return(data)
}  # close no loci removed

if(length(loci_to_filter)>0){

filt_cr<-data[,-which(colnames(data) %in% loci_to_filter)]

if(nrow(md_df)==length(loci_to_filter)) print("all loci were removed")
if(nrow(md_df)==length(loci_to_filter)) return(filt_cr)

if(nrow(md_df)>length(loci_to_filter)){

filt_cr<-tidy.df(filt_cr)
head(filt_cr)

print(paste("Original data: ",nrow(data)," individuals; ",length(data)-2," loci",sep=""))
print(paste("Filtered data: ",nrow(filt_cr)," individuals; ",length(filt_cr)-2," loci",sep=""))
print(paste((length(data)-2)-(length(filt_cr)-2)," loci with more than ", max_missing*100," % missing data removed",sep=""))

return(list(miss_summ=md_df,filt_dat=filt_cr))

} # close not all loci removed

} # close some loci removed

} # close missing_data

allele_freq<-function(data){

# data = a dataframe of dart SNPs, cleaned and formatted in STEP_01

# based on script from evachan.org, modified for dartseq's unusual coding of heterozygotes:
# 0 = Reference allele homozygote,	e.g. AA
# 1 = SNP allele homozygote,		e.g. BB
# 2 = heterozygote					e.g. AB

## n, n0, n1, n2: number of samples with total non-missing genotype, and geno=0,1,or 2
## p: allele frequency
## maf & mgf: minor allele & genotype frequencies

geno<-data[,grep("L",colnames(data))[1]:ncol(data)]
ghead(geno)

m <- ncol(geno)     ## number of snps
n <- nrow(geno)     ## number of individuals

## calc_n (genotype frequencies)
n0 <- apply(geno==0,2,sum,na.rm=T)
n1 <- apply(geno==1,2,sum,na.rm=T)
n2 <- apply(geno==2,2,sum,na.rm=T)
n <- n0 + n1 + n2

## calculate allele frequencies (the key difference for dartseq data is n2 in the eqn for p. The heterozygote is often coded as n1 in other data types and in most basic pop gen explanations)
p <- ((2*n0)+n2)/(2*n)
q <- 1 - p

return(list(p=p,q=q,n0=n0,n1=n1,n2=n2,n=n))

} # close allele_freq

# minor allele frequency:
maf_summary<-function(data){

# data = a dataframe of dart SNPs, cleaned and formatted in STEP_01

afreqs<-allele_freq(data)

# based on script from evachan.org, modified for dartseq's unusual coding of heterozygotes:
maf <- pmin(afreqs$p, afreqs$q)
mgf <- apply(cbind(afreqs$n0,afreqs$n1,afreqs$n2),1,min) / afreqs$n
maf<-data.frame(locus=names(maf),maf=as.numeric(maf))
maf<-tidy.df(maf)
return(maf)

} # close maf_summary

# filter loci below a minor allele frequency cut-off:
maf_filter<-function(maf_sum,data,min_maf){

# maf_sum = the summary data frame created with maf_summary
# data = a dataframe of dart SNPs, cleaned and formatted in STEP_01
# min_maf = the minimum maf allowed (between 0:1) 

loci_to_filter<-as.character(maf_sum$locus[which(maf_sum$maf<min_maf)])
filt_maf<-data[,-which(colnames(data) %in% loci_to_filter)]
filt_maf<-tidy.df(filt_maf)

print(paste("Original data: ",nrow(data)," individuals; ",length(data)-2," loci",sep=""))
print(paste("Filtered data: ",nrow(filt_maf)," individuals; ",length(filt_maf)-2," loci",sep=""))
print(paste((length(data)-2)-(length(filt_maf)-2)," loci with minor allele freq < ",min_maf*100," % removed",sep=""))

return(filt_maf)

} # close maf_filter

# Hardy-Weinberg Equilibrium test:
hwe_exact<-function(data){

# data = a dataframe of dart SNPs, cleaned and formatted in STEP_01

afreqs<-allele_freq(data)

# based on script from evachan.org, modified for dartseq's unusual coding of heterozygotes:

## HWE: Fisher's Exact test
# z <- cbind(n0, ceiling(n1/2), floor(n1/2), n2) # orig; updated for dartseq
z <- cbind(afreqs$n0, ceiling(afreqs$n2/2), floor(afreqs$n2/2), afreqs$n1)
z <- lapply(split(z, 1:nrow(z)), matrix, ncol=2)
z <- lapply(z, fisher.test)
hwe.fisher <- as.numeric(unlist(lapply(z, "[[", "estimate")))
hwe.fisher.p <- as.numeric(unlist(lapply(z, "[[", "p.value")))
hwe.fisher.df<-data.frame(locus=names(afreqs$n0),p=hwe.fisher.p)
hwe.fisher.df$p.adj<-p.adjust(hwe.fisher.df$p,method="hochberg")

return(hwe.fisher.df)

} # close hwe_exact
