
# Author: Annabel Smith

# This is much faster than format_genepop (only a few mins for 70,000 markers and 450 individuals)

format_structure<-function(data,headline){

# data = filtered snp data
# headline = used as file name, something descriptive
# depends = gp_param in format_genepop to write parameter file

# Uses similar framework to genepop converter:
gp_cons1<-data[,3:length(data)]
ghead(gp_cons1)

# Re-code genotypes:
# This only works if the is.na goes in the first place:
gp_cons1<-apply(
gp_cons1,2,function(x)
ifelse(is.na(x),"-9 -9",
ifelse(x=="2","1 2",
ifelse(x=="1","2 2",
ifelse(x=="0","1 1",x
))))
)
ghead(gp_cons1)

# Add site and individual:
gp_cons1<-data.frame(data[,c(2,1)],gp_cons1)
ghead(gp_cons1)

# Make locus index:
str_loci<-colnames(data[,3:length(data)])
head(str_loci)

str.file.name<-paste(headline,".txt",sep="")

# Write loci to first line, space separated:
write(paste(str_loci,collapse=" "),file=str.file.name,append=T)

# Then append data:
write.table(gp_cons1,file=str.file.name,append=T,row.names=F,col.names=F,quote=F,sep=" ")

# Write parameters, from format_genepop library:
gp_param(data,headline)

} # close format_structure

