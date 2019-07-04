
# ------------------------------------ #
# ----------- SUPPLEMENT 03  --------- #
# ------------------------------------ #

### LD TESTS
### Author: Annabel Smith

# Load and tidy workspace and remove everything except necessary objects:
load("../04_workspaces/STEP01_proc_wksp"); rm(list=setdiff(ls(), c("snp_onerow","linf","sdat")))

# load functions:
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Start with full data set, snp_onerow

# Remove cultivars and outgroups:
filt_ogcult<-snp_onerow[-c(grep("OG",snp_onerow$site),grep("CAT",snp_onerow$site),grep("CCT",snp_onerow$site),grep("CTP",snp_onerow$site)),]
filt_ogcult<-tidy.df(filt_ogcult)

# Remove any site with < 8 individuals
filt_sml<-filt_ogcult[-which(filt_ogcult$site %in% names(table(filt_ogcult$site)[table(filt_ogcult$site)<8])),]
filt_sml<-tidy.df(filt_sml)
ghead(filt_sml)
sml_nomono<-mono_loci(filt_sml,3)
ghead(sml_nomono)

# randomly sample three sites and 20 loci for testing:
samp_rows<-which(sml_nomono$site %in% sample(levels(sml_nomono$site),3))
samp_cols<-which(colnames(sml_nomono)[3:ncol(sml_nomono)] %in% sample(colnames(sml_nomono)[3:ncol(sml_nomono)],20))
rand_snp<-sml_nomono[samp_rows,c(1,2,samp_cols)]
rand_snp<-tidy.df(rand_snp)
ghead(rand_snp)

# make sure they all have > 8 individuals:
table(rand_snp$site)

# format genotype file for the calc_LD function from evachang and run test:

sites_to_test<-levels(sml_nomono$site)
dat_test<-sml_nomono
start_site<-which(sites_to_test=="CDF")+1

param_file<-paste("LD_results/parameters","_LD.txt",sep="")

for (i in start_site:length(sites_to_test)){

max_miss<-0.2

site.thisrun<-sites_to_test[i]
data.thisrun<-dat_test[dat_test$site==site.thisrun,]
data.thisrun<-tidy.df(data.thisrun)
orig.noloci<-ncol(data.thisrun[,3:ncol(data.thisrun)])
print(site.thisrun,sep="/n/n")
print(Sys.time())
data.thisrun<-mono_loci(data.thisrun,3)
aftermono.noloci<-ncol(data.thisrun[,3:ncol(data.thisrun)])

data.thisrun<-missing_data(data.thisrun,3,max_miss)
aftermiss.noloci<-ncol(data.thisrun[,3:ncol(data.thisrun)])

if(ncol(data.thisrun)==0) out.list[[i]]<-data.frame(site=site.thisrun,loc1="no_loci",loc2="no_loci",r2=NA) 
if(ncol(data.thisrun)==0) next

data.thisrun<-data.thisrun[,3:ncol(data.thisrun)]

# geno needs to be m x n where m is the number of markers and n is the number of individuals:
data.thisrun<-t(data.thisrun)

ld.thisrun<-calc_LD(data.thisrun,inds=1:nrow(data.thisrun), get.D=F, get.Dprime=F, get.rsq=T, get.chisq=F, get.chisq_prime=F)

loc_combn<-combn(rownames(data.thisrun),2)

df_now<-data.frame(site=site.thisrun,loc1=loc_combn[1,],loc2=loc_combn[2,],r2=ld.thisrun$rsq[lower.tri(ld.thisrun$rsq)])

df_now<-df_now[which(df_now$r2>0.5),]
df_now<-tidy.df(df_now)

if(nrow(df_now)==0) df_now<-data.frame(n1=c("no loci under r2=0.5"))

file_name<-paste("LD_results/",site.thisrun,"_LD.txt",sep="")
write.table(df_now,file_name,quote=F,row.names=F,sep=	"\t")

# write parameters:
write(paste("site = ",site.thisrun,sep=""),file=param_file,append=T)
write(as.character(Sys.time()),file=param_file,append=T)
write(paste("Original data: ",nrow(data.thisrun)," individuals; ",orig.noloci," loci",sep=""),file=param_file,append=T)
write(paste("Filtered data: ",nrow(data.thisrun)," individuals; ",aftermono.noloci," loci",sep=""),file=param_file,append=T)
write(paste(orig.noloci-aftermono.noloci," loci monomorphic removed",sep=""),file=param_file,append=T)

write(paste("Original data: ",nrow(data.thisrun)," individuals; ",aftermono.noloci," loci",sep=""),file=param_file,append=T)
write(paste("Filtered data: ",nrow(data.thisrun)," individuals; ",aftermiss.noloci," loci",sep=""),file=param_file,append=T)
write(paste(aftermono.noloci-aftermiss.noloci," loci with more than ", max_miss*100," % missing data removed",sep=""),file=param_file,append=T)

} # close for

# Create summary file from LD_results:
files.now<-dir("LD_results")[3:length(dir("LD_results"))]
names.now<-substr(files.now,1,unlist(gregexpr(".txt",files.now))-1)

# This takes a long time (~20 mins):
ld_sum<-write.table(data.frame(locus_pair=paste(unlist(lapply(files.now,function(x) as.character(read.table(paste("LD_results/",x,sep=""),header=T)[,2]))),unlist(lapply(files.now,function(x) as.character(read.table(paste("LD_results/",x,sep=""),header=T)[,3]))),sep=""),correlation=unlist(lapply(files.now,function(x) as.character(read.table(paste("LD_results/",x,sep=""),header=T)[,4])))), file="LD_summary.txt", sep="\t", row.names=F, quote=F)

# Summarise file dimensions:

xx<-as.numeric(as.character(unlist(lapply(files.now,function(x) length(read.table(paste("LD_results/",x,sep=""),header=T)[,2])))))

file.dimensions<-data.frame(file=names.now,no_rows=xx)
# write.table(file.dimensions,"file_dimensions.txt",quote=F,row.names=F,sep="\t")

# Analyse LD results:

dir("..")
dir("LD_parameters")

occ<-read.table("LD_parameters/LD_summary.txt",header=T)
head(occ)

# Set LD correlation cut-off
ld_cut<-0.75
occ75<-occ[occ$correlation>ld_cut,]
occ75<-tidy.df(occ75)
dim(occ75)
head(occ75)

# Summarise the number of populations in which each locus pair was in LD:
occ_tab<-table(occ75$locus_pair)
occ_summ<-data.frame(locus_pair=names(occ_tab),occ=as.numeric(occ_tab))
head(occ_summ)

occ_summ3<-occ_summ[occ_summ$occ>2,]
occ_summ3<-tidy.df(occ_summ3)
head(occ_summ3)

hist(occ_summ3$occ[occ_summ3$occ<20])
table(occ_summ3$occ)

# Get locus pairs which were consistently in LD in 6 or more populations:
occ_pop6<-occ_summ[occ_summ$occ>5,]
occ_pop6<-tidy.df(occ_pop6)
head(occ_pop6)

# The determine which loci would need to be removed to ensure no linked loci would occur together:

# Add separate loci:
second_L<-unlist(lapply(gregexpr("L",occ_pop6$locus_pair),function(x)x[2]))
occ_pop6$loc1<-substr(occ_pop6$locus_pair,1,(second_L-1))
occ_pop6$loc2<-substr(occ_pop6$locus_pair,second_L,nchar(as.character(occ_pop6$locus_pair)))
check.rows(occ_pop6)
head(occ_pop6)

# The total number of loci in the linked data set:
all_ldloc<-as.character(c(occ_pop6$loc1,occ_pop6$loc2))
length(unique(all_ldloc))
head(all_ldloc)

freq_loc<-data.frame(locus=names(table(all_ldloc)),no_times_total=as.numeric(table(all_ldloc)))
head(freq_loc)
head(occ_pop6)

removed.loci<-list()

for (i in 1:nrow(occ_pop6)){

line.thisrun<-occ_pop6[i,]
l1.thisrun<-line.thisrun$loc1
l2.thisrun<-line.thisrun$loc2

freq_l1<-freq_loc$no_times_total[which(freq_loc$locus==l1.thisrun)]
freq_l2<-freq_loc$no_times_total[which(freq_loc$locus==l2.thisrun)]

# If one of the pair has already been assigned to the "remove" pile, then the pair is OK and can skip to the next test
if(length(which(c(l1.thisrun,l2.thisrun) %in% unlist(removed.loci)==T))>0) next

# If they have the same frequency in the linkage summary, remove l2:
if (freq_l1==freq_l2) {
	removed.loci[[i]]<-l2.thisrun
	next
	} # close if same freq

# If they have a different frequency, remove the one with the higher frequency:
if(freq_l1!=freq_l2) {
	removed.loci[[i]]<-c(l1.thisrun,l2.thisrun)[which(c(freq_l1,freq_l2)==max(freq_l1,freq_l2))]
	} # close different freq

} # close for

rm_loci<-data.frame(locus=unlist(removed.loci))
rm_loci$for_removal<-1
loci_toremove<-as.character(rm_loci$locus)
head(loci_toremove)

# write.table(rm_loci,"LD_r75_over5pop.txt",sep="\t",row.names=F,quote=F)

# Check:
occ_rm<-merge(occ_pop6,rm_loci,by.x="loc1",by.y="locus",all.x=T, all.y=F)
colnames(occ_rm)[length(colnames(occ_rm))]<-"l1_removal"
occ_rm$l1_removal[which(is.na(occ_rm$l1_removal))]<-0
occ_rm<-merge(occ_rm,rm_loci,by.x="loc2",by.y="locus",all.x=T, all.y=F)
colnames(occ_rm)[length(colnames(occ_rm))]<-"l2_removal"
occ_rm$l2_removal[which(is.na(occ_rm$l2_removal))]<-0

# Sometimes both will be removed, because they both occur somewhere else
check.rows(occ_rm)
head(occ_rm)

# But every line should add to at least 1:
range(rowSums(occ_rm[,which(colnames(occ_rm) %in% c("l1_removal","l2_removal"))]))

# This shows that at least one of each pair will be removed and that no linked pair will be included in the final data set. 

save.image("../04_workspaces/Supp03_LD_tests")










