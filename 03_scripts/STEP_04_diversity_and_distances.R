
# ------------------------------------ #
# ------------- STEP 04  ------------- #
# ------------------------------------ #

### Diversity & distances
### Author: Annabel Smith

# Load workspace:
load("../04_workspaces/STEP04_divdist_wksp")

# load functions:
invisible(lapply(paste("../02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Load libraries:
library("diveRsity")
library("geosphere")
library("hierfstat")
library("adegenet")

#########################################
##     GENIND object & site data:      ##
#########################################
{

gp_dir<-"../ANALYSIS_RESULTS/Genepop_DATA_FILES"
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
head(sdat,3); dim(sdat)

#########################################
##     Genetic & enviro distances:     ##
#########################################
{

### -- *** CALCULATE FST:

gp_dir<-"../ANALYSIS_RESULTS/Genepop_DATA_FILES"
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

## -- ** Admixture Diversity Score:

### ---- level of admixture per population

# Working from equations in Harismendy et al. 2019, Journal of the American Medical Informatics Association, 26(5), 457–46

# Structure summary objects loaded from: Supplement_03_plot_structure.R

# 491 individuals including cultivars and outgroups
head(all_dat); dim(all_dat)
head(site_assig); dim(all_dat)

# Associated site data (61 pops including cultivars and outgroups):
head(sdat,2); dim(sdat)
head(sdat2); dim(sdat)

# The same data but with outgroups and cultivars removed (53 pops, 454 individuals):
head(sdpie,2); dim(sdpie)
head(cluster.data); dim(cluster.data)

# Aggregated data to give the assignment proportion for each population (53 pops):
head(r1); dim(r1)
head(agg_dat); dim(agg_dat)

# Use individual level data for the analysis (6 clusters with assignment probabilities for 454 individuals):
ind6<-cluster.data[,grep("assig",colnames(cluster.data))]
head(ind6); dim(ind6)

# The Diversity Score from Harismendy et al. We don't need eqn 1 (the cumulative admixture fraction), because ours already sum to 1. 

# The diversity score is the same as the shannon diversity index but with a scaling factor to account for the number of clusters (Hmax). That is, it makes the diversity score relative to complete evenness. 

DS<-function(x,K) {
	Hmax<-K*((1/K)*log(1/K))
	(-sum(x*log(x)))/-Hmax
	}

# SET K:
K<-6

# Calculate admixture diversity for all individuals:
ds_div<-apply(ind6, 1, DS, K=6)
range(ds_div)

ind6_dv<-data.frame(indiv=cluster.data$indiv, admix_div=ds_div)
ind_site<-cluster.data[,c("site_code","indiv","native","region","country","latitude","longitude")]
ind6_dv<-merge(ind6_dv, ind_site, by="indiv", all.x=T, all.y=F)
head(ind6_dv); dim(ind6_dv)

# Summarise site-level admixture diversity:
addiv_site<-aggregate(admix_div~site_code,FUN=mean,data=ind6_dv)
ex_dat<-ind6_dv[,c("site_code","native","region","country","latitude","longitude")]
ex_dat<-ex_dat[-which(duplicated(ex_dat$site_code)),]
ex_dat<-tidy.df(ex_dat)
head(ex_dat); dim(ex_dat)
addiv_site<-merge(addiv_site, ex_dat, by="site_code", all.x=T, all.y=F)
head(addiv_site); dim(addiv_site)

# Allelic richness:
ar_dat<-read.table(paste("../01_data/gen_div.txt",sep=""),header=T)

# Check all are true:
addiv_site$site_code %in% ar_dat$site

# Add admixture diversity
ar_dat<-merge(ar_dat, addiv_site, by.x="site", by.y="site_code", all.x=T, all.y=F)
head(ar_dat); dim(ar_dat)

# write.table(ar_dat, file="ardat.txt", quote=F, row.names=F, sep="\t")

# Is there a difference in admixture between the native and non-native ranges?
library(lmerTest)
library(lme4)
library(AICcmodavg)

head(ind6_dv); dim(ind6_dv)

dv_ind<-lmer(admix_div~native+(1|site_code), data=ind6_dv)
summary(dv_ind)
admix_df<-predictSE(dv_ind, newdata=data.frame(native=c("native","non_native")), se.fit=T)
admix_pred<-data.frame(native=c("native","non_native"),fit=admix_df$fit, se=admix_df$se.fit)
admix_pred$lci<-admix_pred$fit-(1.96*admix_pred$se)
admix_pred$uci<-admix_pred$fit+(1.96*admix_pred$se)

### PLOT for SI:
{

quartz("",4,4)
par(mfrow=c(2,2),mar=c(4,4,1,1), mgp=c(2.5,0.8,0))

plot(1:2,admix_pred$fit, xlim=c(0.5,2.5), ylim=c(min(admix_pred$lci),max(admix_pred$uci)), pch=20, xaxt="n", xlab="", ylab="admixture diversity", las=1,cex.axis=0.85,cex.lab=0.85)
arrows(1:2, admix_pred$lci, 1:2, admix_pred$uci, length=0.02, code=3, angle=90)
axis(side=1, at=c(1,2), labels=c("native","non-native"), cex.axis=0.85)
mtext(bquote("range "~italic("p = ")~.(round(summary(dv_ind)$coefficients[2,5],3))), side=3, line=0,cex=0.7, adj=0)
mtext("(a)",side=3, line=0, adj=0, at=-0.6, cex=0.8)

blankplot()

plot(ar_dat$ar[ar_dat$native=="native"], ar_dat$admix_div[ar_dat$native=="native"], pch=20, col="black", cex.axis=0.8, las=1,ylab="admixture diversity",xlab="", cex.lab=0.85)
points(ar_dat$ar[ar_dat$native=="non_native"], ar_dat$admix_div[ar_dat$native=="non_native"], pch=20, col="red", cex.axis=0.8, las=1,cex.lab=0.8)
title(xlab="allelic richness (neutral)", cex.lab=0.85, mgp=c(2,1,0))
mtext("(b)",side=3, line=0, adj=0, at=1.119, cex=0.8)

# Is there are r.ship between ar and admixture diversity?
ar_admix<-lm(admix_div~ar*native,data=ar_dat)
mtext(bquote("range x ar"~italic("p = ")~.(round(summary(ar_admix)$coefficients[4,4],3))), side=3, line=0,cex=0.7, adj=0)

plot(ar_dat$ar_adapt[ar_dat$native=="native"], ar_dat$admix_div[ar_dat$native=="native"], pch=20, col="black", cex.axis=0.8, las=1,ylab="admixture diversity",xlab="", cex.lab=0.8)
points(ar_dat$ar_adapt[ar_dat$native=="non_native"], ar_dat$admix_div[ar_dat$native=="non_native"], pch=20, col="red", cex.axis=0.8, las=1,cex.lab=0.8)
title(xlab="allelic richness (adaptive)", cex.lab=0.85, mgp=c(2,1,0))
mtext("(c)",side=3, line=0, adj=0, at=1.121, cex=0.8)

# Is there are r.ship between ar_adapt and admixture diversity?
adapt_admix<-lm(admix_div~ar_adapt*native,data=ar_dat)
mtext(bquote("range x ar"~italic("p = ")~.(round(summary(adapt_admix)$coefficients[4,4],3))), side=3, line=0,cex=0.7, adj=0)

} # close plot

} # close genetic diversity

























