
# ------------------------------------ #
# ------------- STEP 05  ------------- #
# ------------------------------------ #

### environmental analysis
### Author: Annabel Smith

# load functions:
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

load("../04_workspaces/STEP05_env_wksp")
# library("ecodist")
library("gdm")
library("AICcmodavg")
library("adegenet")
library("poppr")

#########################################
####  	     	  DATA: 	   		 ####
#########################################
{

dat_dir<-"../01_data"
dir(dat_dir)

# **** ---- site data:
sdat<-read.table(paste(dat_dir,"site_data.txt",sep="/"),header=T)

# Remove cultivars, outgroups, RCH and sites with small sample sizes:
sdat<-sdat[-which(sdat$site_code=="RCH"),]
sdat<-sdat[-grep("OG",sdat$site_code),]
sdat<-sdat[-which(sdat$site_code %in% c("CAT","CCT","CTP")),]
sml_samp<-as.character(sdat$site_code[which(sdat$n_gt<7)])
sdat<-sdat[-which(sdat$site_code %in% sml_samp),]
sdat<-tidy.df(sdat)
head(sdat,2); dim(sdat)

# **** ---- pairwise pop stats:
pw<-read.table(paste(dat_dir,"pw_pop_stats.txt",sep="/"),header=T)
head(pw)

# **** ---- REMOVE OUTGROUPS, CULTIVARS and sites with small sample sizes:
cult_lines<-c(which(pw$pop1 %in% c("CAT","CCT","CTP")),which(pw$pop2 %in% c("CAT","CCT","CTP")))
og_lines<-unique(c(grep("OG",pw$pop1),grep("OG",pw$pop2)))
sml_lines<-unique(c(which(pw$pop1 %in% sml_samp),which(pw$pop2 %in% sml_samp)))

pw<-pw[-unique(c(og_lines,cult_lines,sml_lines)),]
pw<-tidy.df(pw)
head(pw); dim(pw)

save.image("../04_workspaces/STEP05_env_wksp")

# **** ---- site level genetic diversity:
gd<-read.table(paste(dat_dir,"gen_div.txt",sep="/"),header=T)
head(gd); dim(gd)

# **** ---- REMOVE OUTGROUPS, CULTIVARS and sites with small sample sizes:
gd<-gd[-grep("OG",gd$site),]
gd<-gd[-which(gd$site %in% c("CAT","CCT","CTP")),]
gd<-gd[-which(gd$max_n<7),]
gd<-tidy.df(gd)
head(gd,2); dim(gd)

# **** ---- indiv level genetic diversity:
ind<-read.table(paste(dat_dir,"indiv_het.txt",sep="/"),header=T)
ind<-ind[-grep("OG",ind$site),]
ind<-ind[-which(ind$site %in% c("CAT","CCT","CTP")),]
ind<-ind[which(ind$site %in% gd$site),]
ind$ind_het<-as.numeric(as.character(ind$ind_het))
ind<-tidy.df(ind)
head(ind)

# Add average individual het to genetic diversity data:
gd$ih_avg<-tapply(ind$ind_het,ind$site,mean)
head(gd,3); dim(gd)

# ** -- Add genetic diversity data site data:
gds<-merge(sdat,gd,by.x="site_code",by.y="site",all.x=T,all.y=F)
gds$max_n<-NULL
head(gds,3); dim(gds)

save.image("../04_workspaces/STEP05_env_wksp")

# **** ---- Add distances to pw pop stats:
head(pw,3)
head(gds,3)

# Site data for sites with genetic data (for GDM):
sdt<-sdat[-which(sdat$genetics=="N"),]
sdt<-tidy.df(sdt)
head(sdt,2); dim(sdt)
head(pw,2)

# These should all be true:
sdt$site_code %in% unique(c(as.character(pw$pop1),as.character(pw$pop2)))

### *** To add back later: ,"Vcov","Bground","Vheight", "elevation"
# GDM doesn't work with NAs and the biotic vars have NAs

sll<-sdt[,c("site_code","mt","st","ap","sp","mm","sm")]

head(pw)
head(sll)
m1df<-merge(pw,sll,by.x="pop1",by.y="site_code",all.x=T,all.y=F)
colnames(m1df)[colnames(m1df) %in% c("mt","st","ap","sp","mm","sm")]<-c("mt1","st1","ap1","sp1","mm1","sm1")
m1df<-merge(m1df,sll,by.x="pop2",by.y="site_code",all.x=T,all.y=F)
colnames(m1df)[colnames(m1df) %in% c("mt","st","ap","sp","mm","sm")]<-c("mt2","st2","ap2","sp2","mm2","sm2")
m1df<-m1df[order(m1df$pop1,m1df$pop2),]
m1df<-tidy.df(m1df)
m1df<-m1df[,c(2,1,3:length(m1df))]
head(m1df,2)

# Add distances:
pw$tmp_dist<-abs(m1df$mt1-m1df$mt2)
pw$tmpv_dist<-abs(m1df$st1-m1df$st2)
pw$precip_dist<-abs(m1df$ap1-m1df$ap2)
pw$precipv_dist<-abs(m1df$sp1-m1df$sp2)
pw$moist_dist<-abs(m1df$mm1-m1df$mm2)
pw$moistv_dist<-abs(m1df$sm1-m1df$sm2)

# Remove unused cols:
pw<-pw[,-which(colnames(pw) %in% c("lat1","lon1","lat2","lon2"))]
head(pw,2)

# Check all altitudes etc. are present
pw$pop1[which(is.na(pw$moist_dist))]

# Make Euro subset:
europop<-sdt$site_code[which(sdt$region=="Europe")]
epw<-pw
epw<-epw[which(epw$pop1 %in% europop),]
epw<-epw[which(epw$pop2 %in% europop),]
head(epw,2)

# Make non-native subset:
nnpop<-sdt$site_code[which(sdt$region!="Europe")]
nnpw<-pw
nnpw<-nnpw[which(nnpw$pop1 %in% nnpop),]
nnpw<-nnpw[which(nnpw$pop2 %in% nnpop),]
head(nnpw,2)

save.image("../04_workspaces/STEP05_env_wksp")

# **** ---- Set modeling variables:

demo_resp<-c("Y0_ros_m2","emp_pgr_ros","ros_reprod_m2")

enviro_pred<-colnames(gds)[colnames(gds) %in% c("mt","ap","sp","st")]

gen_resp<-colnames(gds)[colnames(gds) %in% c("He","ar","ar_adapt")]

all_resp<-c(demo_resp,gen_resp)

gds<-gds[,c(1:which(colnames(gds)=="longitude"),which(colnames(gds) %in% c(demo_resp,enviro_pred,gen_resp)))]
head(gds,2); dim(gds)

# **** ---- Add meta data:
env_rows<-which(!is.na(gds$mt))
dem_rows<-which(!is.na(gds$Y0_ros_m2))
gen_rows<-which(!is.na(gds$He))

gds$dem_env<-ifelse(rownames(gds) %in% as.numeric(names(table(c(env_rows,dem_rows))[which(table(c(env_rows,dem_rows))==2)])),1,0)
gds$gen_dem<-ifelse(rownames(gds) %in% as.numeric(names(table(c(gen_rows,dem_rows))[which(table(c(gen_rows,dem_rows))==2)])),1,0)
gds$gen_env<-ifelse(rownames(gds) %in% as.numeric(names(table(c(gen_rows,env_rows))[which(table(c(gen_rows,env_rows))==2)])),1,0)

gds[,c("site_code","Y0_ros_m2","He","ar","mt","dem_env","gen_dem","gen_env")]

# **** ---- SCALE environmental predictors:
gdsc<-cbind(gds[,1:which(colnames(gds)=="longitude")],apply(gds[,c("mt","ap","sp","st")],2,scale), gds[,which(colnames(gds)=="Y0_ros_m2"):length(gds)])
head(gdsc,3); dim(gdsc)
head(gds,3); dim(gds)

# Check distributions of responses:
{
quartz(title="",width=6,height=4, dpi=100)
par(mfrow=c(2,3),mgp=c(2.5,1,0),oma=c(0,0,0,0),mar=c(5,4,2,1))

for (i in 1:length(all_resp)){
resp.thisrun<-all_resp[i]
hist(gds[,resp.thisrun],main=resp.thisrun,xlab="")
} # close for

# Genetic responses OK
# Emp pgr OK

# NOT OK:
# density
# reprod_m2

# LOG offenders:
to_log<-colnames(gdsc)[c(grep("Y0",colnames(gdsc)),grep("reprod_m2",colnames(gdsc)))]
dont_log<-colnames(gdsc)[-which(colnames(gdsc) %in% to_log)]
gdsc<-cbind(gdsc[,dont_log],apply(gdsc[,to_log],2,function(x) log(x+1)))
head(gdsc,3)

# ID logged
resp.df<-data.frame(resp=all_resp)
resp.df<-data.frame(resp=all_resp,logged=ifelse(all_resp %in% to_log, 1,0))

# Re-check:
quartz(title="",width=6,height=4, dpi=100)
par(mfrow=c(2,3),mgp=c(2.5,1,0),oma=c(0,0,0,0),mar=c(5,4,2,1))

for (i in 1:length(all_resp)){

resp.thisrun<-all_resp[i]

if(resp.df$logged[i]==0) main.thisrun<-resp.thisrun else main.thisrun<-paste("log(",resp.thisrun,"+1)",sep="")

hist(gdsc[,resp.thisrun],main= main.thisrun,xlab="")

} # close for

} # close check response

save.image("../04_workspaces/STEP05_env_wksp")	
} # close data

#########################################
####  	    DATA EXPLANATIONS: 	   	 ####
#########################################

# sdt: environment & genetic data for GDM:
# nrow=53
head(sdt,3); dim(sdt)

# gds: UNSCALED DATA genetic diversity, site, environment and demography:
# nrow=63
head(gds,3); dim(gds)

# gdsc: SCALED DATA genetic diversity, site, environment and demography, non-normal responses LOGGED:
# nrow=63
head(gdsc,3); dim(gdsc)

#########################################
## 				 AMOVA				   ##
#########################################
{
	
###-- genind objects:
# See parameter files in gp_dir for filters
genind_filt1 # no OG or cultivars
genind_filt2 # all sites
genind_filt3 # non-neutral loci, no OG or cultivars, no small sample sizes

###-- Define strata:

# Do this for population genind_filt1@pop and for range (native/non_native):

# First, change pop names to match sdat, so we can add native/non-native

strata(genind_filt1)<- data.frame(pop=pop(genind_filt1))

strata(genind_filt1)$pop<-substr(strata(genind_filt1)$pop,1,nchar(as.character(strata(genind_filt1)$pop))-1)
new.rnames<-as.character(strata(genind_filt1)$pop)

new.rnames[which(new.rnames=="VIR")]<-"VA"

new.rnames[grep("_",substr(new.rnames,nchar(new.rnames),nchar(new.rnames)))]<-substr(new.rnames[grep("_",substr(new.rnames,nchar(new.rnames),nchar(new.rnames)))],1,nchar(new.rnames[grep("_",substr(new.rnames,nchar(new.rnames),nchar(new.rnames)))])-1)

strata(genind_filt1)$pop<-new.rnames
strata(genind_filt1)$pop

# Check that all these are true:
table(as.character(strata(genind_filt1)$pop) %in% sdat$site_code)

site_range_df<-sdat[which(sdat$site_code %in% as.character(strata(genind_filt1)$pop)),c("site_code","native")]
colnames(site_range_df)<-c("pop","range")
head(site_range_df); dim(site_range_df)

strata(genind_filt1)<-merge(strata(genind_filt1),site_range_df,by="pop",all.x=T,all.y=F)
head(strata(genind_filt1)); dim(strata(genind_filt1))

# Need to set dist=NULL as it doesn't run with it set on F 
amova1<-poppr.amova(genind_filt1,~range/pop,dist=NULL, within=T)

amova2<-poppr.amova(genind_filt1,~range/pop,dist=NULL, within=F)

amova3<-poppr.amova(genind_filt1,~pop,dist=NULL, within=F)

amova1 # within individuals
amova2 # within populations
amova3 # within ranges

save.image("../04_workspaces/Supp_amova")

} # close amova

#########################################
####  	     	  GDM: 	 	  		 ####
#########################################
{

# Check for correlations among envrionmental variables:
{

# Full data set:
{
head(sdt,3); dim(sdt)

env<-sdt[,c("native","mt","st","ap","sp","mm","sm")]
# make native/non_native == 0/1
env$native<-ifelse(env$native=="native",0,1)
head(env)

# Multi-collinearity
mcl1<-mcl_v3("env")
mcl1<-mcl1[order(-mcl1$VIF),]
mcl1

# Remove until all < 3:
env2<-env[,-which(colnames(env)=="mm")]
head(env2)
mcl2<-mcl_v3("env2")
mcl2<-mcl2[order(-mcl2$VIF),]
mcl2

env3<-env2[,-which(colnames(env2)=="sm")]
head(env3)
mcl3<-mcl_v3("env3")
mcl3<-mcl3[order(-mcl3$VIF),]
mcl3

# Pair-wise correlations:
env.cor<-cor(env,use="na.or.complete")
env.names<-combn(colnames(env),2)

quartz("",10,6,dpi=100,pointsize=8)
# mgp here controls ylab; xlab is controlled in the function
par(mfrow=c(5,6),mar=c(4,4,2,0),oma=c(0,0,1,2),mgp=c(2.8,1,0))
plot.pw("env",env.cor,env.names)

} # close all correlations

# EUR only:
{

seur<-sdat[sdat$region=="Europe",]
seur<-tidy.df(seur)
head(seur,3)

env_eur<-seur[,c("mt","st","ap","sp","mm","sm")]
head(env_eur)

# Multi-collinearity
eur_mcl1<-mcl_v3("env_eur")
eur_mcl1<-eur_mcl1[order(-eur_mcl1$VIF),]
eur_mcl1

# Remove until all < 3:
eur_env2<-env_eur[,-which(colnames(env_eur)=="mm")]
head(eur_env2)
eur_mcl2<-mcl_v3("eur_env2")
eur_mcl2<-eur_mcl2[order(-eur_mcl2$VIF),]
eur_mcl2

# sm is not a problem for EUR, but removing it for consistency with all other analyses in the paper:
eur_env3<-eur_env2[,-which(colnames(eur_env2)=="sm")]
head(eur_env3)
eur_mcl3<-mcl_v3("eur_env3")
eur_mcl3<-eur_mcl3[order(-eur_mcl3$VIF),]
eur_mcl3

# Pair-wise correlations:
eur_env.cor<-cor(env_eur,use="na.or.complete")
eur_env.names<-combn(colnames(env_eur),2)

quartz("",10,6,dpi=100,pointsize=8)
# mgp here controls ylab; xlab is controlled in the function
par(mfrow=c(4,6),mar=c(4,4,2,0),oma=c(0,0,1,2),mgp=c(2.8,1,0))
plot.pw("env_eur",eur_env.cor,eur_env.names)

} # close euro correlations

# NN only:
{

snn<-sdat[sdat$region!="Europe",]
snn<-tidy.df(snn)
head(snn,3)

env_nn<-snn[,c("mt","st","ap","sp","mm","sm")]
head(env_nn)

# Multi-collinearity
nn_mcl1<-mcl_v3("env_nn")
nn_mcl1<-nn_mcl1[order(-nn_mcl1$VIF),]
nn_mcl1

# Remove until all < 3:
nn_env2<-env_nn[,-which(colnames(env_nn)=="sm")]
head(nn_env2)
nn_mcl2<-mcl_v3("nn_env2")
nn_mcl2<-nn_mcl2[order(-nn_mcl2$VIF),]
nn_mcl2

nn_env3<-nn_env2[,-which(colnames(nn_env2)=="mm")]
head(nn_env3)
nn_mcl3<-mcl_v3("nn_env3")
nn_mcl3<-nn_mcl3[order(-nn_mcl3$VIF),]
nn_mcl3

# Pair-wise correlations:
nn_env.cor<-cor(env_nn,use="na.or.complete")
nn_env.names<-combn(colnames(env_nn),2)

quartz("",10,6,dpi=100,pointsize=8)
# mgp here controls ylab; xlab is controlled in the function
par(mfrow=c(4,6),mar=c(4,4,2,0),oma=c(0,0,1,2),mgp=c(2.8,1,0))
plot.pw("env_nn",nn_env.cor,nn_env.names)

} # close nn correlations

} # close correlations

### FULL DATA SET:
{

# Individual enviro variables:

# The data format for GDM is close to the merged data (m1df) that was used above to create the pw_pop_stats. 
head(m1df,2)
# Remove "nat_dist" in the colname:
m2df<-m1df
m2df<-m2df[,-grep("nat_dist",colnames(m2df))]
head(m2df,2)
mcl3 # variable correlations

# NOTES: This data set includes all variables except three that were removed for multi-collinearity or correlations. Removing the worst VIF offenders - sm and mm - was sufficient to reduce VIF to < 3 (in all data sets: all pops, eur and nn). 
vars.toincludeALL<-as.character(mcl3$response)
vars.toincludeALL[vars.toincludeALL=="native"]<-"nat"

# Input the name of the first column containing enviro vars, i.e. not including the fsts and geog coords:
enviro.indALL<-which(colnames(m2df)=="nat1")

# Check that there are only two of each in the colnames of the relevant data set. This should return a 2 x nvar matrix:
sapply(vars.toincludeALL,function(x)grep(x,colnames(m2df)[enviro.indALL:length(colnames(m2df))]))

gdm_all<-format_gdm(m2df,vars.toincludeALL,enviro.indALL)
head(gdm_all,3)

} # close full data set

### EUROPE ONLY:
{

# Individual enviro variables:

# Add EUR index to m1:
head(m1df,2)
head(sdat,2)
sr<-sdat[,c("site_code","region")]

meur<-merge(m1df,sr,by.x="pop1",by.y="site_code")
colnames(meur)[which(colnames(meur)=="region")]<-"region1"
meur<-merge(meur,sr,by.x="pop2",by.y="site_code")
colnames(meur)[which(colnames(meur)=="region")]<-"region2"
meur<-meur[,c(2,1,3:length(meur))]

meur<-meur[meur$region1=="Europe" & meur$region2=="Europe",]
meur<-tidy.df(meur)
head(meur,3)

# Remove all Switzerland except the lowest elevation:

swiss.lines<-as.character(sdat[grep("SW",sdat$site_code),"site_code"][-which(sdat[grep("SW",sdat$site_code),"elevation"]==min(sdat[grep("SW",sdat$site_code),"elevation"]))])

meur<-meur[-unique(c(which(meur$pop1 %in% swiss.lines),which(meur$pop2 %in% swiss.lines))),]
meur<-tidy.df(meur)

# Remove mm to reduce VIF < 3. Mean temperature can stay in. 
vars.toincludeEUR<-as.character(eur_mcl3$response)
eur_mcl3 # variable correlations

# Input the name of the first column containing enviro vars, i.e. not including the fsts and geog coords:
enviro.indEUR<-which(colnames(meur)=="mt1")

# Check that there are only two of each in the colnames of the relevant data set. This should return a 2 x nvar matrix:
sapply(vars.toincludeEUR,function(x)grep(x,colnames(meur)[enviro.indEUR:length(colnames(meur))]))

gdm_eur<-format_gdm(meur,vars.toincludeEUR,enviro.indEUR)
head(gdm_eur,3)

} # close EUR

### NON-NATIVE ONLY:
{

# Individual enviro variables:

# Add NN index to m1:
mnn<-merge(m1df,sr,by.x="pop1",by.y="site_code")
colnames(mnn)[which(colnames(mnn)=="region")]<-"region1"
mnn<-merge(mnn,sr,by.x="pop2",by.y="site_code")
colnames(mnn)[which(colnames(mnn)=="region")]<-"region2"
mnn<-mnn[,c(2,1,3:length(mnn))]

mnn<-mnn[mnn$region1!="Europe" & mnn$region2!="Europe",]
mnn<-tidy.df(mnn)
head(mnn,3)

# For the NN data set, the only variable that needs to be removed is mean moisture to reduce VIF < 3. Mean temperature can stay in. But remove BG as it's highly correlated with veg cover. 
vars.toincludeNN<-as.character(nn_mcl3$response)
nn_mcl3 # variable correlations

# Input the name of the first column containing enviro vars, i.e. not including the fsts and geog coords:
enviro.indNN<-which(colnames(mnn)=="mt1")

# Check that there are only two of each in the colnames of the relevant data set. This should return a 2 x nvar matrix:
sapply(vars.toincludeNN,function(x)grep(x,colnames(mnn)[enviro.indNN:length(colnames(mnn))]))

gdm_nn<-format_gdm(mnn,vars.toincludeNN,enviro.indNN)
head(gdm_nn,3)

} # close non-native

save.image("../04_workspaces/STEP05_env_wksp")

# Fit GDMs (indiv enviro vars) -----------------------------------------
{

# Fit:
{
GEO=T # use geographic distance as a predictor?

head(gdm_all,3)
head(gdm_eur,3)
head(gdm_nn,3)

class(gdm_all)<-c("gdmData","data.frame")
class(gdm_eur)<-c("gdmData","data.frame")
class(gdm_nn)<-c("gdmData","data.frame")

# Fit models
vars.toincludeALL # all sites, all vars
gdmALL <- gdm(data=gdm_all, geo=GEO)
summary(gdmALL)
str(gdmALL)
gdmALL$coefficients
gdmALL$explained

gdmEUR <- gdm(gdm_eur, geo=GEO)
summary(gdmEUR)
str(gdmEUR)
gdmEUR$coefficients
gdmEUR$explained

gdmNN <- gdm(gdm_nn, geo=GEO)
summary(gdmNN)
str(gdmNN)
gdmNN$coefficients
gdmNN$explained

# extract spline data for custom plotting
allSplines <- isplineExtract(gdmALL)
eurSplines <- isplineExtract(gdmEUR)
nnSplines <- isplineExtract(gdmNN)

allSplines$x[1:5,] # 7 variables
allSplines$y[1:5,]

eurSplines$x[1:5,] # 6 variables
eurSplines$y[1:5,]

nnSplines$x[1:5,] # 6 variables
nnSplines$y[1:5,]

save.image("../04_workspaces/STEP05_env_wksp")

} # close fit

# xlabs:
{

gdm_xlabs_all<-data.frame(spline=c("Geographic","nat","mt","st","ap","sp"),newname=c("Geographic distance","Origin (native/non-native)","Mean temperature","Temperature seasonality","Mean precipitation (mm)","Precipitation seasonality (mm)"))

all_xlabs<-format_gdm_xlabs(allSplines)
eur_xlabs<-format_gdm_xlabs(eurSplines)
nn_xlabs<-format_gdm_xlabs(nnSplines)

} # close xlabs

# Test importance
{

allTest <- gdm.varImp(gdm_all, geo=T, nPerm=50, parallel=T, cores=10)
eurTest <- gdm.varImp(gdm_eur, geo=T, nPerm=50, parallel=T, cores=10)
nnTest <- gdm.varImp(gdm_nn, geo=T, nPerm=50, parallel=T, cores=10)

save.image("../04_workspaces/STEP05_env_wksp")

# Deviance explained and model p-values:
allTest[[1]] 
eurTest[[1]]
nnTest[[1]]

# Individual predictor p-values:
allTest[[3]][,1]
eurTest[[3]][,1]
nnTest[[3]][,1]

p.adjust(allTest[[3]][,1],method="bonferroni")
p.adjust(eurTest[[3]][,1],method="bonferroni")
p.adjust(nnTest[[3]][,1],method="bonferroni")

} # close test

} # close fit models

#### PLOT ------------------------------------------

# Don't use plot_gdm, there is an unfixed bug, use plot_gdm_simp, or do it manually

# Individual environmental vars:
{
### FULL DATA SET:

quartz("",8,6,dpi=90)
par(mfrow=c(3,4),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm(allSplines,gdmALL,all_xlabs,allTest)

# simplified for SI:
quartz("",6,4,dpi=90)
par(mfrow=c(2,3),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm_simp(allSplines,gdmALL,all_xlabs,allTest,adjust_p=T)

### EUROPE & NON-NATIVE:

quartz("",8,6,dpi=90)
par(mfrow=c(3,4),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm(eurSplines,gdmEUR,eur_xlabs,eurTest)

# simplified for SI:
quartz("",6,4,dpi=90)
par(mfrow=c(2,3),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm_simp(eurSplines,gdmEUR,eur_xlabs,eurTest,adjust_p=T)

quartz("",8,6,dpi=90)
par(mfrow=c(3,4),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm(nnSplines,gdmNN,nn_xlabs,nnTest)

# simplified for SI:
quartz("",6,4,dpi=90)
par(mfrow=c(2,3),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm_simp(nnSplines,gdmNN,nn_xlabs,nnTest,adjust_p=T)

} # close plot indiv vars

# Plot for paper:
{

# Plotting raw FSTs in main doc, splines separately in SI

# euro data:
head(eurSplines$y); dim(eurSplines$y)
head(meur,2); dim(meur)
head(gdm_eur,2); dim(gdm_eur)
head(eurSplines$x)
head(eurSplines$y)
eurTest[[3]][,1][which(rownames(eurTest[[3]])=="Geographic")]

# Set spline variables for plotting:
x.now<-eurSplines$x[,"Geographic"]/10
y.now<-eurSplines$y[,"Geographic"]

x.sp<-eurSplines$x[,"sp"]
y.sp<-eurSplines$y[,"sp"]
sp_dist<-abs(meur$sp1-meur$sp2)

quartz("",8,4,dpi=80,pointsize=18)
par(mfrow=c(1,2),mar=c(3,4,1,0),oma=c(0,0,1,1),mgp=c(2,0.8,0),xpd=NA)

# (a) Europe geographic
plot(x.now,y.now,ylim=c(range(meur$fst)[1],0.5),cex=0.5,pch=20,xlab="",ylab=expression("Genetic distance ("*italic("F")[ST]*")"),type="n",las=1,col="red",lwd=2,cex.axis=0.7,cex.lab=0.8, xlim=c(0,3),xaxt="n")

title(xlab="Geographic distance (,000 km)",mgp=c(1.6,1,0),cex.lab=0.8)
mtext("(a)",side=2, las=1,line=1, at=0.57, cex=0.8, adj=0)
points((meur$geog_dist/100000)/10, meur$fst,pch=20, cex=0.3)
axis(side=1, at=0:3, labels=0:3,cex.axis=0.7, mgp=c(2.3,0.5,0))

# (b) Europe precip
plot(x.sp,y.sp,ylim=c(range(meur$fst)[1],0.5),cex=0.5,pch=20,xlab="",ylab=expression("Genetic distance ("*italic("F")[ST]*")"),type="n",las=1,col="red",lwd=2,xlim=range(min(c(sp_dist)),max(c(sp_dist))),cex.axis=0.7,cex.lab=0.8,xaxt="n")

title(xlab="Precipitation seasonality (mm)",mgp=c(1.6,1,0),cex.lab=0.8)
mtext("(b)",side=2, las=1,line=1, at=0.57, cex=0.8, adj=0)
points(sp_dist,meur$fst,pch=20, cex=0.3)
axis(side=1, cex.axis=0.7, mgp=c(2.3,0.5,0))

} # close plot for paper

} # close GDM

#########################################
##   Gen div ~ environ + demography:   ##
#########################################
{

### DATA SET UP:
{

# MAIN DATA:
head(gds,2); dim(gds)
head(gdsc,2); dim(gds)

# Use same base data for all three stages:

# Stage 1: demog ~ environment
# Stage 1: genetic ~ environment
# Stage 2: genetic ~ demography

# nobs==42 (without outliers)
# emp_pgr is NOT comparable with AIC (fewer obs), but we're only using AIC for within responses, not across models

# unscaled:
st1<-gds[which(gds$gen_dem==1),]
st1<-tidy.df(st1)

# scaled:
st1sc<-gdsc[which(gdsc$gen_dem==1),]
st1sc<-tidy.df(st1sc)

# UNLOG demographic predictors and SCALE:
demo_resp
st1sc[,c("Y0_ros_m2","ros_reprod_m2")]<-st1[,c("Y0_ros_m2","ros_reprod_m2")]
st1sc[,which(colnames(st1sc) %in% demo_resp)]<-apply(st1sc[,which(colnames(st1sc) %in% demo_resp)],2,scale)

# Check for outliers in predictors:
preds_all<-c(enviro_pred,demo_resp)

quartz(title="",width=6,height=4, dpi=90, pointsize=14)
par(mfrow=c(2,4),mgp=c(2.5,1,0),oma=c(0,0,0.5,0),mar=c(2,4,1,1.5))

for (i in 1:length(preds_all)){
pred.now<-preds_all[i]
boxplot(st1sc[,pred.now],main=pred.now)
} # close for

# Remove single reprod outlier:
reprod.outl<-which(st1sc$ros_reprod_m2>3)
st1sc<-st1sc[-reprod.outl,]
st1sc<-tidy.df(st1sc)
st1<-st1[-reprod.outl,]
st1<-tidy.df(st1)

# Remove single density outlier:
dens.outl<-which(st1sc$Y0_ros_m2>2)
st1sc<-st1sc[-dens.outl,]
st1sc<-tidy.df(st1sc)
st1<-st1[-dens.outl,]
st1<-tidy.df(st1)

# gdsc has scaled environmental predictors and unscaled demographic responses. Non-normal responses are logged: use for demog ~ environment:
env_sc<-gdsc[which(gdsc$site_code %in% st1sc$site_code),]
env_sc<-tidy.df(env_sc)
env_unsc<-gds[which(gds$site_code %in% st1$site_code),]
env_unsc<-tidy.df(env_unsc)

} # close data set-up

### MODEL SET UP:

# Stage 1: demog ~ environment
# Stage 1: genetic ~ environment 
# Stage 2: genetic ~ demography

# DATA FOR GENETIC RESPONSES:
head(st1,2); dim(st1) # nothing scaled
head(st1sc,2); dim(st1sc) # all predictors scaled, use for genetic ~ environment and genetic ~ demography

# DATA FOR DEMOGRAPHIC RESPONSES:
# environmental predictors scaled, demographic responses not scaled, use for demog ~ environment:
head(env_sc,3); dim(env_sc) 
head(env_unsc,3); dim(env_unsc)

demo_resp # demographic response
gen_resp # genetic response

preds_all # all predictors
enviro_pred # enviro predictors

mod_tab<-data.frame(type=c(rep("dem_env",12),rep("gd_env",8)),response=c(rep(demo_resp,each=4),rep(gen_resp[2:3],each=4)),predictor=enviro_pred,sc_data=c(rep("env_sc",12),rep("st1sc",8)), unsc_data=c(rep("env_unsc",12),rep("st1",8)))
mod_tab<-apply(mod_tab, 2, as.character)

mod_tab2<-data.frame(type=rep("gd_dem",6),response=rep(gen_resp[2:3],each=3), predictor=demo_resp,sc_data=rep("st1sc",6),unsc_dat=rep("st1",6))
mod_tab2<-apply(mod_tab2, 2, as.character)
mod_tab2

mod_tab<-data.frame(rbind(mod_tab,mod_tab2))
mod_tab$mod_name<-paste(mod_tab$response, mod_tab$predictor, sep="_")
head(mod_tab)

# Scaled data
head(env_sc,3); dim(env_sc) 
head(st1sc,2); dim(st1sc)

### RUN MODELS:

aic.store<-list()
coef.store<-list()

# Run models, rank by AICc, store AIC results and top model coefficients:

for (i in 1:nrow(mod_tab)){

resp.thisrun<-as.character(mod_tab$response[i])
pred.thisrun<-as.character(mod_tab$predictor[i])
data.thisrun<-get(as.character(mod_tab$sc_data[i]))
head(data.thisrun,2); dim(data.thisrun)

mod_null<-lm(get(resp.thisrun)~1,data=data.thisrun)
mod_nat<-lm(get(resp.thisrun)~native,data=data.thisrun)
mod_pred<-lm(get(resp.thisrun)~get(pred.thisrun),data=env_sc)
mod_add<-lm(get(resp.thisrun)~get(pred.thisrun)+native,data=env_sc)
mod_int<-lm(get(resp.thisrun)~get(pred.thisrun)*native,data=env_sc)

# AICc:
# If delta is negative, it didn't improve the model, for positive and negative AIC values
aic.df<-data.frame(response=resp.thisrun,predictor=pred.thisrun,mod_name=paste(resp.thisrun, pred.thisrun,sep="_"),model=c("null","nat","pred","add","int"),AICc=c(AICc(mod_null), AICc(mod_nat), AICc(mod_pred), AICc(mod_add), AICc(mod_int)),rank="unsupported")
aic.df$rank<-as.character(aic.df$rank)

aic.df$rank[which(aic.df$AICc==min(aic.df$AICc))]<-"best"

best.model<-as.character(aic.df$model[aic.df$rank=="best"])

aic.ofnull<-aic.df$AICc[which(aic.df$model=="null")]

aic.ofbest<-aic.df$AICc[which(aic.df$rank=="best")]
aic.df$delta.ofbest<-aic.ofbest-aic.df$AICc

# does best beat null?
# if not, then null becomes best:

if(best.model!="null" & aic.ofnull-aic.ofbest<2){

aic.df$rank[which(aic.df$rank=="best")]<-"unsupported"
aic.df$rank[which(aic.df$model=="null")]<-"best"
best.model<-as.character(aic.df$model[aic.df$rank=="best"])

}

if(best.model=="null" & length(which(aic.df$delta.ofbest>=2))>0) aic.df$rank[which(aic.df$delta.ofbest>=2)]<-"supported"

within.best<-NULL
best.index<-NULL
within.best<-which(abs(aic.df$delta.ofbest)<=2)
best.index<-which(aic.df$rank=="best")
within.best<-within.best[-which(within.best==best.index)]

if(best.model!="null" & length(within.best)>0) {

aic.df$rank[within.best]<-"supported"
	
	}

aic.store[[i]]<-aic.df

# store coefficients for all supported models:

supp.mods<-unique(as.character(aic.df$model[-which(aic.df$rank=="unsupported")]))

coef.temp<-list()

for (j in 1:length(supp.mods)){

m.now<-supp.mods[j]

mod.name<-paste("mod_",m.now,sep="")

coef.now<-summary(get(mod.name))$coefficients
coef.df<-data.frame(Response=resp.thisrun,Predictor=pred.thisrun,rownames(coef.now),coef.now)
colnames(coef.df)<-c("Response","Predictor","Term","Estimate","SE","t_value","p_value")
coef.df<-tidy.df(coef.df)

coef.temp[[j]]<-coef.df

} # close j

coef.tempres<-do.call(rbind,coef.temp)

coef.store[[i]]<-coef.tempres

} # close i

aic.res<-do.call(rbind,aic.store)
head(aic.res)

coef.res<-do.call(rbind,coef.store)
head(coef.res)

# write.table(coef.res, "coef.res.txt",row.names=F, quote=F, sep="\t")

# -------------------- #
# BEST models:
# -------------------- #

### GET ESTIMATES:

best.mods<-aic.res[which(aic.res$rank=="best"),]
best.mods<-best.mods[-which(best.mods$model=="null"),]
best.mods<-best.mods[-which(best.mods$model=="nat"),]
best.mods<-tidy.df(best.mods)
best.mods

mt2<-mod_tab[,c("mod_name","sc_data","unsc_data")]
head(mt2)

best.mods<-merge(best.mods, mt2, by="mod_name")

# arrange for plotting:
best.mods<-best.mods[c(4,3,1,2),]
best.mods<-tidy.df(best.mods)
best.mods

# Get estimates from model for plots:

preds.store<-list()

for (i in 1:nrow(best.mods)){

resp.thisrun<-as.character(best.mods$response[i])
pred.thisrun<-as.character(best.mods$predictor[i])
scdat.thisrun<-get(as.character(best.mods$sc_data[i]))
head(scdat.thisrun,2); dim(scdat.thisrun)
range.thisrun<-NULL

mod.thisrun<-NULL

type.thisrun<-as.character(best.mods$model[i])

if(type.thisrun=="add") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)+native,data=scdat.thisrun)
if(type.thisrun=="int") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)*native,data=scdat.thisrun)

summary(mod.thisrun)

nd.thisrun<-data.frame(native=c(rep("native",50),rep("non_native",50)),pred=seq(min(scdat.thisrun[,pred.thisrun],na.rm=T),max(scdat.thisrun[,pred.thisrun],na.rm=T),length.out=50))
colnames(nd.thisrun)[2]<-pred.thisrun
head(nd.thisrun)
pr1<-predict(mod.thisrun,newdata=nd.thisrun,se.fit=T)
pr2<-data.frame(nd.thisrun,fit=pr1$fit,se=pr1$se.fit)
pr2$lci<-pr2$fit-(1.96*pr2$se)
pr2$uci<-pr2$fit+(1.96*pr2$se)
colnames(pr2)[2]<-"predictor"

dfpred<-data.frame(pred=pred.thisrun,pr2)
head(dfpred)

preds.store[[i]]<-dfpred

if(is.null(range.thisrun)==F){

summary(mod.thisrun)

nd2<-data.frame(native=c(rep("native",1),rep("non_native",1)),pred=mean(scdat.thisrun[,pred.thisrun],na.rm=T))
colnames(nd2)[2]<-pred.thisrun
pr3<-predict(mod.thisrun,newdata=nd2,se.fit=T)
pr4<-data.frame(nd2,fit=pr3$fit,se=pr3$se.fit)
pr4$lci<-pr4$fit-(1.96*pr4$se)
pr4$uci<-pr4$fit+(1.96*pr4$se)
colnames(pr4)[2]<-"predictor"
pr4<-data.frame(resp=resp.thisrun,pred=pred.thisrun,pr4)
head(pr4)

range.store[[i]]<-pr4

} # close range

} # close i

save.image("../04_workspaces/STEP05_env_wksp")

### PLOT:

best.mods # models to plot (n=4)
preds.store # environmental estimates (n=3)

plot.data<-best.mods
plot.data$resp_lab<-c(rep("Population growth rate",1),rep("Neutral genetic diversity",1),rep("Adaptive genetic diversity",2))
plot.data$pred_lab<-c(rep("Mean temperature",1),rep("Temperature seasonality",1),rep("Mean precipitation (mm)",1),rep("Temperature seasonality",1))
plot.data$space_after<-c(rep("no",4))
plot.data$lett<-c("a","b","c","d")

quartz(title="",width=8,height=6.5, dpi=90, pointsize=16.5)
par(mfrow=c(2,2),mgp=c(2.7,0.9,0),oma=c(0,0,0.5,6.5),mar=c(4.5,4,1.5,1.7))

for (i in 1:nrow(plot.data)){

resp.thisrun<-as.character(plot.data$response[i])
pred.thisrun<-as.character(plot.data$predictor[i])
est.thisrun<-preds.store[[i]]

scdat.thisrun<-get(as.character(plot.data$sc_data[i]))
unscdat.thisrun<-get(as.character(plot.data$unsc_data[i]))
head(scdat.thisrun,2)

raw.scaled<-scdat.thisrun[,c(which(colnames(scdat.thisrun) %in% c("site","native","mt","st","ap","sp")),which(colnames(scdat.thisrun)==resp.thisrun))]
raw.unscaled<-unscdat.thisrun[,c(which(colnames(unscdat.thisrun) %in% c("site","native","mt","st","ap","sp")),which(colnames(unscdat.thisrun)==resp.thisrun))]
head(raw.scaled); dim(raw.scaled)
head(raw.unscaled); dim(raw.unscaled)

yl.now<-c(min(c(min(est.thisrun$lci),min(raw.scaled[,resp.thisrun],na.rm=T))),max(c(max(est.thisrun$uci),max(raw.scaled[,resp.thisrun],na.rm=T))))

head(est.thisrun)

plot(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],ylim=yl.now,type="n",bty="l",xlab="",ylab=plot.data$resp_lab[i],xaxt="n", las=1,cex.axis=0.8)

pg.ci("predictor","est.thisrun",x.subset="native",colour=rgb(0,0,0,0.1))
lines(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],lty=1)
lines(est.thisrun[est.thisrun$native=="non_native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="non_native"],lty=2,col="red")

points(raw.scaled[raw.scaled$native=="native",pred.thisrun],raw.scaled[raw.scaled$native=="native",resp.thisrun],col="black",cex=0.5, pch=20)
points(raw.scaled[raw.scaled$native=="non_native",pred.thisrun],raw.scaled[raw.scaled$native=="non_native",resp.thisrun],col="red",cex=0.5, pch=20)

if(resp.thisrun=="emp_pgr_ros") arrows(-10,0,2.5,0,length=0,lwd=0.3,col="black")

xax<-seq(round(min(raw.scaled[,pred.thisrun]),0),round(max(raw.scaled[,pred.thisrun]),0),length.out=5)
head(raw.unscaled)
xaxl<-round(seq(round(min(raw.unscaled[,pred.thisrun]),0),round(max(raw.unscaled[,pred.thisrun]),0),length.out=5),0)

axis(side=1, at=xax, labels=xaxl, cex.axis=0.8)

xlab.thirun<-as.character(plot.data$pred_lab[i])
if(length(grep("emperatu",xlab.thirun))==0) title(xlab= xlab.thirun,mgp=c(2.2,1,0)) else title(xlab=bquote(.(xlab.thirun)*" ("*degree*"C)"),mgp=c(2.2,1,0))

mtext(paste("(",plot.data$lett[i],")",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.9,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/7))

# if(plot.data$space_after[i]=="yes") blankplot()

# Add general legend:
par(xpd=NA)
if(i==2) legend(par("usr")[2]+((par("usr")[2]/10)*1),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("native","non native"),lty=c(1,2),bty="n",col=c("black","red"))
if(i==2) legend(par("usr")[2]+((par("usr")[2]/10)*2.75),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("",""),pch=c(20),col=c("black","red"),bty="n")
par(xpd=F)

} # close plots i

# -------------------- #
# SUPPORTED models:
# -------------------- #

head(aic.res)

supp.mods<-aic.res[which(aic.res$rank=="supported"),]
supp.mods<-supp.mods[-which(supp.mods$model=="nat"),]
supp.mods<-tidy.df(supp.mods)
supp.mods

supp.mods<-merge(supp.mods, mt2, by="mod_name")

# arrange for plotting:
supp.mods<-supp.mods[c(7,6,5,3,4,2,1),]
supp.mods<-tidy.df(supp.mods)
supp.mods

# Get estimates from model for plots:

preds.store.supp<-list()

for (i in 1:nrow(supp.mods)){

resp.thisrun<-as.character(supp.mods$response[i])
pred.thisrun<-as.character(supp.mods$predictor[i])
scdat.thisrun<-get(as.character(supp.mods$sc_data[i]))
head(scdat.thisrun,2); dim(scdat.thisrun)
range.thisrun<-NULL

mod.thisrun<-NULL

type.thisrun<-as.character(supp.mods$model[i])

if(type.thisrun=="add") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)+native,data=scdat.thisrun)
if(type.thisrun=="int") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)*native,data=scdat.thisrun)

summary(mod.thisrun)

nd.thisrun<-data.frame(native=c(rep("native",50),rep("non_native",50)),pred=seq(min(scdat.thisrun[,pred.thisrun],na.rm=T),max(scdat.thisrun[,pred.thisrun],na.rm=T),length.out=50))
colnames(nd.thisrun)[2]<-pred.thisrun
head(nd.thisrun)
pr1<-predict(mod.thisrun,newdata=nd.thisrun,se.fit=T)
pr2<-data.frame(nd.thisrun,fit=pr1$fit,se=pr1$se.fit)
pr2$lci<-pr2$fit-(1.96*pr2$se)
pr2$uci<-pr2$fit+(1.96*pr2$se)
colnames(pr2)[2]<-"predictor"

dfpred<-data.frame(pred=pred.thisrun,pr2)
head(dfpred)

preds.store.supp[[i]]<-dfpred

if(is.null(range.thisrun)==F){

summary(mod.thisrun)

nd2<-data.frame(native=c(rep("native",1),rep("non_native",1)),pred=mean(scdat.thisrun[,pred.thisrun],na.rm=T))
colnames(nd2)[2]<-pred.thisrun
pr3<-predict(mod.thisrun,newdata=nd2,se.fit=T)
pr4<-data.frame(nd2,fit=pr3$fit,se=pr3$se.fit)
pr4$lci<-pr4$fit-(1.96*pr4$se)
pr4$uci<-pr4$fit+(1.96*pr4$se)
colnames(pr4)[2]<-"predictor"
pr4<-data.frame(resp=resp.thisrun,pred=pred.thisrun,pr4)
head(pr4)

range.store[[i]]<-pr4

} # close range

} # close i

save.image("../04_workspaces/STEP05_env_wksp")

### PLOT:

supp.mods # models to plot (n=7)
preds.store.supp 

pd.supp<-supp.mods
pd.supp$resp_lab<-c(rep("Population growth rate",1),rep("Neutral genetic diversity",4),rep("Adaptive genetic diversity",2))
pd.supp$pred_lab<-c(rep("Mean temperature",1),rep("Density",1),rep("Temperature seasonality",1),rep("Mean precipitation (mm)",2),rep("Mean temperature",1),rep("Mean precipitation (mm)",1))
pd.supp$space_after<-c(rep("no",1),"yes",rep("no",5))
pd.supp$lett<-letters[1:7]

quartz(title="",width=9.5,height=9.5, dpi=70, pointsize=21)
par(mfrow=c(3,3),mgp=c(2.7,0.9,0),oma=c(0,0,0.5,1),mar=c(4.5,4,1.5,1.7))

for (i in 1:nrow(pd.supp)){

resp.thisrun<-as.character(pd.supp$response[i])
pred.thisrun<-as.character(pd.supp$predictor[i])
est.thisrun<-preds.store.supp[[i]]
head(est.thisrun)

scdat.thisrun<-get(as.character(pd.supp$sc_data[i]))
unscdat.thisrun<-get(as.character(pd.supp$unsc_data[i]))
head(scdat.thisrun,2)

raw.scaled<-scdat.thisrun[,c(which(colnames(scdat.thisrun) %in% c("site","native","mt","st","ap","sp","Y0_ros_m2")),which(colnames(scdat.thisrun)==resp.thisrun))]
raw.unscaled<-unscdat.thisrun[,c(which(colnames(unscdat.thisrun) %in% c("site","native","mt","st","ap","sp","Y0_ros_m2")),which(colnames(unscdat.thisrun)==resp.thisrun))]
head(raw.scaled); dim(raw.scaled)
head(raw.unscaled); dim(raw.unscaled)

yl.now<-c(min(c(min(est.thisrun$lci),min(raw.scaled[,resp.thisrun],na.rm=T))),max(c(max(est.thisrun$uci),max(raw.scaled[,resp.thisrun],na.rm=T))))

head(est.thisrun)

plot(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],ylim=yl.now,type="n",bty="l",xlab="",ylab=pd.supp$resp_lab[i],xaxt="n", las=1,cex.axis=0.8)

pg.ci("predictor","est.thisrun",x.subset="native",colour=rgb(0,0,0,0.1))
lines(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],lty=1)
lines(est.thisrun[est.thisrun$native=="non_native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="non_native"],lty=2,col="red")

head(raw.scaled)
points(raw.scaled[raw.scaled$native=="native",pred.thisrun],raw.scaled[raw.scaled$native=="native",resp.thisrun],col="black",cex=0.5, pch=20)
points(raw.scaled[raw.scaled$native=="non_native",pred.thisrun],raw.scaled[raw.scaled$native=="non_native",resp.thisrun],col="red",cex=0.5, pch=20)

# if(resp.thisrun=="emp_pgr_ros") arrows(-10,0,2.5,0,length=0,lwd=0.3,col="black")

xax<-seq(round(min(raw.scaled[,pred.thisrun]),0),round(max(raw.scaled[,pred.thisrun]),0),length.out=5)
head(raw.unscaled)
xaxl<-round(seq(round(min(raw.unscaled[,pred.thisrun]),0),round(max(raw.unscaled[,pred.thisrun]),0),length.out=5),0)

axis(side=1, at=xax, labels=xaxl, cex.axis=0.72)

xlab.thirun<-as.character(pd.supp$pred_lab[i])
if(length(grep("emperatu",xlab.thirun))==0) title(xlab= xlab.thirun,mgp=c(2.2,1,0)) else title(xlab=bquote(.(xlab.thirun)*" ("*degree*"C)"),mgp=c(2.2,1,0))

mtext(paste("(",pd.supp$lett[i],")",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.7,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/7))

# Add general legend:
par(xpd=NA)
if(i==2) legend(par("usr")[2]+((par("usr")[2]/100)*0),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("native","non native"),lty=c(1,2),bty="n",col=c("black","red"))
if(i==2) legend(par("usr")[2]+((par("usr")[2]/100)*20),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("",""),pch=c(20),col=c("black","red"),bty="n")
par(xpd=F)

if(pd.supp$space_after[i]=="yes") blankplot()

} # close plots i


# -------------------- #
# CORRELATION plots:
# -------------------- #

# mt, st, ap

env.testdf<-env_unsc[,c("mt","st","ap")]
head(env.testdf)

env.cortest<-cor(env.testdf,use="na.or.complete")
env.cornames<-combn(colnames(env.testdf),2)

axis.names<-data.frame(cornames=unique(paste(env.cornames)),newnames=c("Mean temperature","Temperature seasonality","Mean precipitation (mm)"))

quartz(title="",width=6.5,height=6.5, dpi=90, pointsize=16)
par(mfrow=c(2,2),mgp=c(2.7,0.9,0),oma=c(0,0,0.5,1),mar=c(4.5,4,1.5,1.7),xpd=F)

for (i in 1:length(env.testdf)){

v1<-env.cornames[1,i]
v2<-env.cornames[2,i]

col1<-which(colnames(env_unsc)==v1)
col2<-which(colnames(env_unsc)==v2)

v1<-as.character(axis.names$newnames[which(axis.names$cornames==v1)])
v2<-as.character(axis.names$newnames[which(axis.names$cornames==v2)])

plot(env_unsc[env_unsc$native=="native",col1], env_unsc[env_unsc$native=="native", col2],type="p",bty="o",xlab="",ylab="", las=1,cex.axis=0.9, col="black",cex=0.5, pch=20,xlim=c(min(env_unsc[,col1]),max(env_unsc[,col1])), ylim=c(min(env_unsc[, col2]),max(env_unsc[, col2])))
points(env_unsc[env_unsc$native=="non_native",col1], env_unsc[env_unsc$native=="non_native",col2],col="red",pch=20, cex=0.5)

if(length(grep("emperatu",v1))==0) title(xlab=v1,mgp=c(2.2,1,0)) else title(xlab=bquote(.(v1)*" ("*degree*"C)"),mgp=c(2.2,1,0))
if(length(grep("emperatu",v2))==0) title(ylab=v2,mgp=c(2.8,1,0)) else title(ylab=bquote(.(v2)*" ("*degree*"C)"),mgp=c(2.2,1,0))

mtext(paste("(",letters[i],")",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.9,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/10))

st_mt<-cor.test(env_unsc[,col1], env_unsc[,col2])
text(par("usr")[2]-((par("usr")[2]/100)*4),(par("usr")[4]/6)*5.5,bquote(italic("r")*" = "*.(round(st_mt$estimate,2))),cex=0.7,adj=1)
text(par("usr")[2]-((par("usr")[2]/100)*4),(par("usr")[4]/6)*5,bquote(italic("p")*" = "*.(round(st_mt$p.value,2))),cex=0.7,adj=1)

# Add general legend:
par(xpd=NA)
# if(i==1) legend(par("usr")[2]+((par("usr")[2]/100)*0),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("native","non native"),lty=c(1,2),bty="n",col=c("black","red"))
if(i==1) legend(par("usr")[2]+((par("usr")[2]/100)*20),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("native","non native"),pch=c(20),col=c("black","red"),bty="n")
par(xpd=F)

if(i==1) blankplot()

} # close i

} # close gen div



























