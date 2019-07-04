
# Author: Annabel Smith

# Warning: format_genepop is very SLOW (up to 1.5 hr for 70,000 markers and 450 individuals)

format_genepop<-function(data,headline){

# data = filtered snp data
# headline = used as file name, something descriptive

gp_cons1<-data[,3:length(data)]
ghead(gp_cons1)

# Re-code genotypes:
# This only works if the is.na goes in the first place:
gp_cons1<-apply(
gp_cons1,2,function(x)
ifelse(is.na(x),"0000",
ifelse(x=="2","0102",
ifelse(x=="1","0202",
ifelse(x=="0","0101",x
))))
)
ghead(gp_cons1)

# Add site and individual:
gp_cons1<-data.frame(data[,1:2],gp_cons1)
ghead(gp_cons1)

## ~~~~ ****** Genepop file ****** ~~~~ ##

# *** WARNING: SLOW
# Hours for full data set
# populations must be in the first column called "site" in data.toconvert
# individuals must be in the second column called "ind" in data.toconvert

ghead(gp_cons1)
gp.head<-headline
pg.file.name<-paste(headline,".gen",sep="")
gp.loci<-paste(colnames(gp_cons1)[grep("L",colnames(gp_cons1))],collapse=", ")
sites<-as.character(unique(gp_cons1$site))
no.sites<-length(sites)
line.breaks<-which(!duplicated(gp_cons1$site))

gp_cons1$ind<-paste(gp_cons1$ind,",",sep="")
ghead(gp_cons1)

for (i in 1:no.sites){

pop.thisrun<-gp_cons1[gp_cons1$site==sites[i],]
ghead(pop.thisrun)

if(i==1) write.table(paste(gp.head,gp.loci,sep="\n"),file=pg.file.name,append=F,row.names=F,col.names=F,quote=F)

write("pop",file=pg.file.name,append=T)

write.table(pop.thisrun[,2:length(pop.thisrun)],file=pg.file.name,append=T,row.names=F,col.names=F,quote=F)

} # close for site

gp_param(data,headline)

# write locus info index:
write.table(data.frame(lind=1:length(colnames(data)[3:ncol(data)]),locus=colnames(data)[3:ncol(data)]),paste(headline,"_loci.txt",sep=""),row.names=F,quote=F,sep="\t")

} # close format_genepop function

gp_param<-function(data,headline)
{
cat(headline,
paste("Time = ",Sys.time(),sep=""),
paste("Number of sites = ",param.nosites,sep=""),
paste("Sites included = ",paste(param.sites,collapse=", "),sep=""),
paste("Number of loci = ",param.noloci,sep=""),
paste("Number of individuals = ",param.noindiv,sep=""),
paste("Monomorphic loci removed = ",param.mono,sep=""),
if(param.repavg==T) paste("Minimum reproducibility = ",ra*100,"%",sep="") else paste("No reproducibility filter",sep=""),
if(param.callrate==T) paste("Maximum missing data = ",cr*100,"%",sep="") else paste("No missing data filter",sep=""),
if(param.MAF==T) paste("Minor allele frequency filter = ", malim*100,"%",sep="") else paste("No MAF filter",sep=""),
if(param.LD==T) paste("Linkage disequilibrium filter = ", ldf*100,"%",sep="") else paste("No LD filter",sep=""),
if(param.HWE==T) paste("HWE filter applied = ",param.HWE,sep="") else paste("No HWE filter",sep=""),
if(param.neu==T) paste("Non-neutral markers removed = ",param.neu,sep="") else paste("Non-neutral markers not removed",sep=""),
if(param.dup==T) paste("Duplicate sequences removed = ",param.HWE,sep="") else paste("Duplicate sequences not removed",sep=""),
file=paste(headline,"_param.txt",sep=""),sep="\n",append=T)

} # close gp_param


