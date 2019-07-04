
# ------------------------------------ #
# ------------- STEP 04  ------------- #
# ------------------------------------ #

### Diversity & distances
### Author: Annabel Smith

# UPDATE 22 March 2019: all updated after bioclim update

# Load and tidy workspace and remove everything except necessary objects:
load("../04_workspaces/STEP01_proc_wksp"); rm(list=setdiff(ls(), c("linf","sdat")))

# Load working workspace:
load("../04_workspaces/STEP04_divdist_wksp")

# load functions:
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Load libraries:
library("adegenet")
library("diveRsity")
library("geosphere")
library("hierfstat")

#########################################
##     GENIND object & site data:      ##
#########################################
{

gp_dir<-"/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSIS_RESULTS/Genepop_DATA_FILES"
dir(gp_dir)

# Make genind objects:
# ~~
genind_filt1<-read.genepop(file=paste(gp_dir,"genepop_filt1.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_filt1
save.image("../04_workspaces/STEP04_divdist_wksp")

# ~~
genind_filt2<-read.genepop(file=paste(gp_dir,"genepop_filt2.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_filt2
save.image("../04_workspaces/STEP04_divdist_wksp")

# ~~ this is the non-neutral data set:
genind_filt3<-read.genepop(file=paste(gp_dir,"genepop_filt3.gen",sep="/"), ncode=2L,quiet=FALSE)
genind_filt3
save.image("../04_workspaces/STEP04_divdist_wksp")

} # close genind

# RESULT:
# See parameter files in gp_dir for filters
genind_filt1 # no OG or cultivars
genind_filt2 # all sites
genind_filt3 # non-neutral loci, no OG or cultivars, no small sample sizes
head(sdat,3)

#########################################
##     Genetic & enviro distances:     ##
#########################################
{

### -- *** CALCULATE FST:

gp_dir<-"/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSIS_RESULTS/Genepop_DATA_FILES"
dir(gp_dir)

# USE ALL POPULATIONS, including outgroups and cultivars (genepop_filt2.gen); can subset this later, but we need all the FSTs:

# Get FST:
# 1hr 20min for 513 x 18320
print(Sys.time())
fst<-diffCalc(paste(gp_dir,"genepop_filt2.gen",sep="/"),fst=T,pairwise=T)
print(Sys.time())
save.image("../04_workspaces/STEP04_divdist_wksp")

head(fst$pairwise$Fst)
head(fst$pairwise$gst)
head(fst$pairwise$Gst)
head(fst$pairwise$GGst)
head(fst$pairwise$D[,1:5])

fstn<-rownames(fst$pairwise$Fst)
fstn<-substr(fstn,1,(nchar(fstn)-2))

ind_endnum<-which(unlist(gregexpr("[0-9]",substr(fstn,nchar(fstn),nchar(fstn))))>0)
with_endnum<-fstn[ind_endnum]
fstn[ind_endnum]<-substr(with_endnum,1,nchar(with_endnum)-1)

ind_endund<-which(unlist(gregexpr("_",substr(fstn,nchar(fstn),nchar(fstn))))>0)
with_endund<-fstn[ind_endund]
fstn[ind_endund]<-substr(with_endund,1,nchar(with_endund)-1)

fst_df<-data.frame(pop1=combn(fstn,2)[1,],pop2=combn(fstn,2)[2,],fst=fst$pairwise$Fst[lower.tri(fst$pairwise$Fst)],gst=fst$pairwise$gst[lower.tri(fst$pairwise$gst)],Gst=fst$pairwise$Gst[lower.tri(fst$pairwise$Gst)],GGst=fst$pairwise$GGst[lower.tri(fst$pairwise$GGst)],D=fst$pairwise$D[lower.tri(fst$pairwise$D)])
head(fst_df)

og_lines<-c(grep("OG",fst_df$pop1),grep("OG",fst_df$pop2))
nonog_lines<-which(rownames(fst_df) %in% og_lines==F)

mean(fst_df$fst[og_lines])
mean(fst_df$fst[nonog_lines])

range(fst_df$fst[og_lines],na.rm=T)
range (fst_df$fst[nonog_lines],na.rm=T)

# write.table(fst_df,"fst.txt",row.names=F,quote=F,sep="\t")

### -- *** ADD GEOGRAPHIC DISTANCE:

dat_dir<-"../01_data"
dir(dat_dir)

# load fst data:
pwpop<-read.table(paste(dat_dir,"pw_pop_stats.txt",sep="/"),header=T)
head(pwpop)

# Update site names:
pwpop$pop1<-as.character(pwpop$pop1)
pwpop$pop2<-as.character(pwpop$pop2)
pwpop$pop1[which(pwpop$pop1=="VIR")]<-"VA"
pwpop$pop2[which(pwpop$pop2=="VIR")]<-"VA"
pwpop$pop1<-as.factor(pwpop$pop1)
pwpop$pop2<-as.factor(pwpop$pop2)

# site data:
head(sdat,3)

# Check all sites in data have site data:
unique(c(levels(pwpop$pop1),levels(pwpop$pop2))) %in% sdat$site_code

sll<-sdat[,c("site_code","latitude","longitude")]
head(sll)
m1<-merge(pwpop,sll,by.x="pop1",by.y="site_code",all.x=T,all.y=F)
colnames(m1)[colnames(m1) %in% c("latitude","longitude")]<-c("lat1","lon1")
m1<-merge(m1,sll,by.x="pop2",by.y="site_code",all.x=T,all.y=F)
colnames(m1)[colnames(m1) %in% c("latitude","longitude")]<-c("lat2","lon2")
m1<-m1[order(m1$pop1,m1$pop2),]
m1<-tidy.df(m1)
m1<-m1[,c(2,1,3:length(m1))]
m1$geog_dist<-distGeo(m1[,c("lon1","lat1")],m1[,c("lon2","lat2")])
m1$lat_dist<-abs(m1$lat1)-abs(m1$lat2)
head(m1)
head(sdat,4)

# Add native distance:
ndis<-sdat[,c("site_code","native")]
head(ndis)
m2<-merge(m1,ndis,by.x="pop1",by.y="site_code",all.x=T,all.y=F)
colnames(m2)[colnames(m2) %in% c("native")]<-c("nat1")
m2<-merge(m2,ndis,by.x="pop2",by.y="site_code",all.x=T,all.y=F)
colnames(m2)[colnames(m2) %in% c("native")]<-c("nat2")
m2<-m2[order(m2$pop1,m2$pop2),]
m2<-tidy.df(m2)
m2<-m2[,c(2,1,3:length(m2))]
head(m2)

m2$nat_dist<-m2$nat1==m2$nat2
m2$nat_dist<-ifelse(m2$nat_dist==T,0,1)
check.rows(m2)

# write.table(m2,"m2.txt",row.names=F,quote=F,sep="\t")

} # close distances

#########################################
##  	  Genetic diversity:  	       ##
#########################################
{

# Use hierfstat functions. divBasic from diveRsity does not work on laptop or big mac with full data set (after hours of running the program quits and says "your computer has run out of application memory"), although it works on test data sets. 
# het<-divBasic("Genepop_files/Development/irl_test")

## -- ** POPULATION LEVEL GENETIC DIVERSITY:

# Remember to update VIR to VA

# Calculate genetic diversity per population in hierfstat (conservative dataset):
gendiv_filt2 <- basic.stats(genind_filt2, diploid = TRUE, digits = 2)
str(gendiv_filt2)
head(gendiv_filt2$Ho)
tail(gendiv_filt2$Ho)
save.image("../04_workspaces/STEP04_divdist_wksp")

gd_filt2<-data.frame(site=names(apply(gendiv_filt2$Ho,2,mean,na.rm=T)),max_n=apply(gendiv_filt2$n.ind.samp,2,max,na.rm=T),Ho=apply(gendiv_filt2$Ho,2,mean,na.rm=T),He=apply(gendiv_filt2$Hs,2,mean,na.rm=T),Fis=apply(gendiv_filt2$Fis,2,mean,na.rm=T))
gd_filt2<-tidy.df(gd_filt2)
head(gd_filt2)

gd_filt2$site<-substr(gd_filt2$site,1,nchar(as.character(gd_filt2$site))-1)

gd_filt2$site[grep("_",substr(gd_filt2$site,nchar(gd_filt2$site),nchar(gd_filt2$site)))]<-substr(gd_filt2$site[grep("_",substr(gd_filt2$site,nchar(gd_filt2$site),nchar(gd_filt2$site)))],1,nchar(gd_filt2$site[grep("_",substr(gd_filt2$site,nchar(gd_filt2$site),nchar(gd_filt2$site)))])-1)

# write.table(gd_filt2,"gd_filt2.txt",sep="\t",row.names=F,quote=F)

## -- ** ALLELIC RICHNESS:

# Remember to update VIR to VA

head(sdat,3); dim(sdat)
range(sdat$n_gt)
sdat[,c("site_code","n_gt")]

# Re-do allelic richness on 53 sites, all with 7-9 individuals:
# See Sept 2018 old code file for old calcs
# Approx. 6 min for 454 x 17162
genind_filt1
genind_filt3

print(Sys.time())
ar_default_rd<- allelic.richness(genind_filt1, diploid = TRUE)
print(Sys.time())
save.image("../04_workspaces/STEP04_divdist_wksp")

print(Sys.time())
ar_default_adapt<- allelic.richness(genind_filt3, diploid = TRUE)
print(Sys.time())
save.image("../04_workspaces/STEP04_divdist_wksp")

# Summarise

# Neutral:
head(ar_default_rd$Ar,2)
ar_default_rd$min.all

ar_res_rd<-data.frame(
site=levels(genind_filt1@pop),
ar_default_rd=apply(ar_default_rd$Ar,2,mean,na.rm=T))

# Non-neutral:
head(ar_default_adapt$Ar,2)
ar_default_adapt$min.all

ar_res_adapt<-data.frame(
site=levels(genind_filt3@pop),
ar_default_adapt=apply(ar_default_adapt$Ar,2,mean,na.rm=T))

test_df<-cbind(ar_res_rd,ar_adapt=ar_res_adapt[,2])
head(test_df)
plot(test_df$ar_default_rd, test_df$ar_adapt)

# write.table(ar_res_adapt,"ar_res_adapt.txt",sep="\t",row.names=F,quote=F)

save.image("../04_workspaces/STEP04_divdist_wksp")

## -- ** INDIVIDUAL HETEROZYGOSITY:

# Calculate individual heterozygosity
# Use dartseq format files, stored in the Genepop folder:
dir(gp_dir)

ds_filt2<-read.table(paste(gp_dir,"dartseq_filt2.txt",sep="/"),header=T)
ghead(ds_filt2); dim(ds_filt2)

# Codes for onerow formatted data:
# 0 = Reference allele homozygote
# 1 = SNP allele homozygote
# 2 = heterozygote

# 3 MINS on laptop

ih_out<-data.frame(ds_filt2[,1:2],ind_het=NA)
head(ih_out,25)

for(i in 1:nrow(ih_out)){

ind.cons<-as.character(ds_filt2[i,3:length(ds_filt2)])
t.cons<-data.frame(ind.cons=as.numeric(names(table(ind.cons))),count=as.numeric(table(ind.cons)),stringsAsFactors=F)

if(length(which(is.na(t.cons$ind.cons)))>0) t.cons<-t.cons[-which(is.na(t.cons$ind.cons)),] else stop("no zero category")

if(length(which(t.cons$ind.cons==2))==0) ih_out[i,3]<-"no_heterozygotes" else ih_out[i,3]<-t.cons$count[t.cons$ind.cons==2]/sum(t.cons$count)

} # close for i
save.image("../04_workspaces/STEP04_divdist_wksp")
head(ih_out)

# write.table(ih_out,"ih.txt",quote=F,row.names=F,sep="\t")

} # close genetic diversity































