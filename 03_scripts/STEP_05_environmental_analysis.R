
# ------------------------------------ #
# ------------- STEP 05  ------------- #
# ------------------------------------ #

### environmental analysis
### Author: Annabel Smith

# load functions:
invisible(lapply(paste("../02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

load("../04_workspaces/STEP05_env_wksp")
library("gdm")
library("adegenet")
library("poppr")
library("AICcmodavg")
library(boot)

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

gen_resp<-colnames(gds)[colnames(gds) %in% c("ar","ar_adapt")]

all_resp<-c(demo_resp,gen_resp)

gds<-gds[,c(1:which(colnames(gds)=="longitude"),which(colnames(gds) %in% c(demo_resp,enviro_pred,gen_resp)))]
head(gds,2); dim(gds)

# **** ---- Add meta data:
env_rows<-which(!is.na(gds$mt))
dem_rows<-which(!is.na(gds$Y0_ros_m2))
gen_rows<-which(!is.na(gds$ar))

gds$dem_env<-ifelse(rownames(gds) %in% as.numeric(names(table(c(env_rows,dem_rows))[which(table(c(env_rows,dem_rows))==2)])),1,0)
gds$gen_dem<-ifelse(rownames(gds) %in% as.numeric(names(table(c(gen_rows,dem_rows))[which(table(c(gen_rows,dem_rows))==2)])),1,0)
gds$gen_env<-ifelse(rownames(gds) %in% as.numeric(names(table(c(gen_rows,env_rows))[which(table(c(gen_rows,env_rows))==2)])),1,0)

head(gds,2); dim(gds)

save.image("../04_workspaces/STEP05_env_wksp")	

} # close data

#########################################
####  	    DATA EXPLANATIONS: 	   	 ####
#########################################

# sdt: environment & genetic data for GDM:
# nrow=53
head(sdt,3); dim(sdt)

# gds: genetic diversity, site, environment and demography, all sites, unscaled:
# nrow=63
head(gds,3); dim(gds)

#########################################
## 				 AMOVA				   ##
#########################################
{

load("../04_workspaces/Supp_amova")
	
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
head(m1df,2); dim(m1df)
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

# bootstrap deviance explained
{

head(gdm_all,3)
head(gdm_eur,3)
head(gdm_nn,3)

# dev explained from the original fit:
gdmALL$explained
gdmEUR$explained
gdmNN$explained

# BOOTSTRAP function for GDM:
stat_gdm <- function (dat, w) {
gdm_test1 <- gdm(dat[w,], geo=GEO)
gdm_test1$explained
}

gdm_all_boot <- boot(gdm_all, stat_gdm, R=10000)
# plot(gdm_all_boot)
boot.ci(gdm_all_boot,type=c("norm","basic","perc","bca"))

gdm_eur_boot <- boot(gdm_eur, stat_gdm, R=10000)
# plot(gdm_eur_boot)
boot.ci(gdm_eur_boot,type=c("norm","basic","perc","bca"))

gdm_nn_boot <- boot(gdm_nn, stat_gdm, R=10000)
# plot(gdm_nn_boot)
boot.ci(gdm_nn_boot, type=c("norm","basic","perc","bca"))

save.image("../04_workspaces/STEP05_env_wksp")

} # close bootstrap

# xlabs:
{

gdm_xlabs_all<-data.frame(spline=c("Geographic","nat","mt","st","ap","sp"),newname=c("Geographic distance","Origin (native/non-native)","Mean temperature","Temperature seasonality","Mean precipitation (mm)","Precipitation seasonality (mm)"))

all_xlabs<-format_gdm_xlabs(allSplines)
eur_xlabs<-format_gdm_xlabs(eurSplines)
nn_xlabs<-format_gdm_xlabs(nnSplines)

} # close xlabs

# Test importance
{

allTest <- gdm.varImp(gdm_all, geo=T, nPerm=500, parallel=T, cores=10)

save.image("../04_workspaces/STEP05_env_wksp")

eurTest <- gdm.varImp(gdm_eur, geo=T, nPerm=500, parallel=T, cores=10)

save.image("../04_workspaces/STEP05_env_wksp")

nnTest <- gdm.varImp(gdm_nn, geo=T, nPerm=500, parallel=T, cores=10)

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

# Individual environmental vars:
{
### FULL DATA SET:

# simplified for SI:
quartz("",6,4,dpi=90)
par(mfrow=c(2,3),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm_simp(allSplines,gdmALL,all_xlabs,allTest,adjust_p=T)

### EUROPE & NON-NATIVE:

# simplified for SI:
quartz("",6,4,dpi=90)
par(mfrow=c(2,3),mar=c(4,5,1,0),oma=c(0,0,0,1),mgp=c(2.3,1,0),xpd=NA)
plot_gdm_simp(eurSplines,gdmEUR,eur_xlabs,eurTest,adjust_p=T)

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

title(xlab="Geographic distance (1,000 km)",mgp=c(1.6,1,0),cex.lab=0.8)
mtext("(a)",side=2, las=1,line=1, at=0.57, cex=0.8, adj=0)
points((meur$geog_dist/100000)/10, meur$fst,pch=20, cex=0.3)
axis(side=1, at=0:3, labels=0:3,cex.axis=0.7, mgp=c(2.3,0.5,0))

# (b) Europe precip
plot(x.sp,y.sp,ylim=c(range(meur$fst)[1],0.5),cex=0.5,pch=20,xlab="",ylab=expression("Genetic distance ("*italic("F")[ST]*")"),type="n",las=1,col="red",lwd=2,xlim=range(min(c(sp_dist)),max(c(sp_dist))),cex.axis=0.7,cex.lab=0.8,xaxt="n")

title(xlab="Precipitation seasonality",mgp=c(1.6,1,0),cex.lab=0.8)
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
head(gds,2); dim(gds)

# Use same base data for all three stages:

# Stage 1: demog ~ environment
# Stage 1: genetic ~ environment
# Stage 2: genetic ~ demography

gds_upd<-gds[which(gds$gen_dem==1),]
gds_upd<-tidy.df(gds_upd)
head(gds_upd,3); dim(gds_upd)

# Log fecundity:
gds_upd$LOG_ros_reprod_m2<-log(gds_upd$ros_reprod_m2+1)

# Scale (normalise) environmental predictors:
gds_upd$mt_sc<-as.numeric(scale(gds_upd$mt))
gds_upd$st_sc<-as.numeric(scale(gds_upd$st))
gds_upd$ap_sc<-as.numeric(scale(gds_upd$ap))
gds_upd$sp_sc<-as.numeric(scale(gds_upd$sp))

# and demographic predictors:
gds_upd$Y0_ros_m2_sc<-as.numeric(scale(gds_upd$Y0_ros_m2))
gds_upd$emp_pgr_ros_sc<-as.numeric(scale(gds_upd$emp_pgr_ros))
gds_upd$LOG_ros_reprod_m2_sc<-as.numeric(scale(gds_upd$LOG_ros_reprod_m2))

head(gds_upd,3); dim(gds_upd)

# suffix "_sc" = scaled predictors:
# gds_upd includes scaled  environmental predictors (and unscaled enviro vars for plotting), unscaled demographic responses (Y0_ros_m2, LOG_ros_reprod_m2, emp_pgr_ros) and scaled demographic predictors (Y0_ros_m2_sc, LOG_ros_reprod_m2_sc, emp_pgr_ros_sc). 

### MODEL SET UP:

# Stage 1: demog ~ environment
# Stage 1: genetic ~ environment 
# Stage 2: genetic ~ demography

head(gds_upd,3); dim(gds_upd)

# demographic response
demo_resp<-c("Y0_ros_m2", "emp_pgr_ros", "LOG_ros_reprod_m2")
# genetic response
gen_resp<-c("ar","ar_adapt")

# enviro predictors
enviro_pred<-paste(c("mt","st","ap","sp"),"_sc",sep="")

# demo predictors
demo_pred<-c("Y0_ros_m2_sc", "emp_pgr_ros_sc","LOG_ros_reprod_m2_sc")
# all predictors
preds_all<-c(enviro_pred, demo_pred)

mod_tab<-data.frame(type=c(rep("dem_env",12),rep("gd_env",8)),response=c(rep(demo_resp,each=4),rep(gen_resp,each=4)),predictor=enviro_pred)

mod_tab<-apply(mod_tab, 2, as.character)

mod_tab2<-data.frame(type=rep("gd_dem",6),response=rep(gen_resp,each=3), predictor=demo_pred)
mod_tab2<-apply(mod_tab2, 2, as.character)
mod_tab2

mod_tab<-data.frame(rbind(mod_tab,mod_tab2))
mod_tab$mod_name<-paste(mod_tab$response, mod_tab$predictor, sep="_")
head(mod_tab)

head(gds_upd,3); dim(gds_upd)

### RUN MODELS:

aic.store<-list()
coef.store<-list()

# Run models, rank by AICc, store AIC results and top model coefficients:

for (i in 1:nrow(mod_tab)){

resp.thisrun<-as.character(mod_tab$response[i])
pred.thisrun<-as.character(mod_tab$predictor[i])
data.thisrun<-gds_upd
head(data.thisrun,2); dim(data.thisrun)

mod_null<-lm(get(resp.thisrun)~1,data=data.thisrun)
mod_nat<-lm(get(resp.thisrun)~native,data=data.thisrun)
mod_pred<-lm(get(resp.thisrun)~get(pred.thisrun),data=data.thisrun)
mod_add<-lm(get(resp.thisrun)~get(pred.thisrun)+native,data=data.thisrun)
mod_int<-lm(get(resp.thisrun)~get(pred.thisrun)*native,data=data.thisrun)

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

aic_wt<-data.frame(aictab(list(mod_null,mod_nat,mod_pred,mod_add,mod_int),modnames=c("null","nat","pred","add","int")))[,c("Modnames","AICcWt","LL","ModelLik","K")]

aic.df<-merge(aic.df, aic_wt, by.x="model", by.y="Modnames", all.x=T, all.y=F)

aic.df<-aic.df[match(c("null","nat","pred","add","int"),aic.df$model),]
aic.df<-tidy.df(aic.df)

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

# write.table(coef.res,"coef.res.txt",row.names=F, quote=F, sep="\t")

# -------------------- #
# BEST models:
# -------------------- #

best.mods<-aic.res[which(aic.res$rank=="best"),]
best.mods<-best.mods[-which(best.mods$model=="null"),]
best.mods<-best.mods[-which(best.mods$model=="nat"),]
best.mods<-tidy.df(best.mods)
best.mods

# Remove the precip seasonality model - we cannot interpret this interaction because the native and non-native ranges don't have comparable values. Include in the SI instead. 
sp_reprod<-best.mods[which(best.mods$predictor=="sp_sc"),]
best.mods<-best.mods[-which(best.mods$predictor=="sp_sc"),]
best.mods<-tidy.df(best.mods)
best.mods

### BOOTSTRAP interaction terms:

library(boot)

# The best models for adaptive genetic diversity were interactive (and for neutral~ap). Would this be the case if native and native datasets had equal sample sizes?

int.mods<-best.mods[which(best.mods$model=="int"),]
int.mods<-tidy.df(int.mods)
int.mods

head(gds_upd,3); dim(gds_upd)

boot_dat1<-gds_upd[,c("site_code","native","st","st_sc","ap","ap_sc","ar","ar_adapt")]
head(boot_dat1,3); dim(boot_dat1)
table(boot_dat1$native)

# BOOTSTRAP function for genetic diversity models:
# Statistic is interactive model coefficient:
stat <- function (dat, w, resp, pred) {
form1<-as.formula(paste(resp,"~",pred,"*native",sep=""))
m1_int<-lm(form1,data=dat[w,])
summary(m1_int)$coefficients[4,1]
}

# Models to bootstrap:
int.mods

ar_ap <- boot(boot_dat1, stat, resp="ar", pred="ap_sc", R=10000, strata=boot_dat1$native)
# plot(ar_ap)
boot.ci(ar_ap, type=c("norm","basic","perc","bca"))

adapt_st <- boot(boot_dat1, stat, resp="ar_adapt", pred="st_sc", R=10000, strata=boot_dat1$native)
# plot(adapt_st)
boot.ci(adapt_st, type=c("norm","basic","perc","bca"))

adapt_ap <- boot(boot_dat1, stat, resp="ar_adapt", pred="ap_sc", R=10000, strata=boot_dat1$native)
# plot(adapt_ap)
boot.ci(adapt_ap, type=c("norm","basic","perc","bca"))
save.image("../04_workspaces/STEP05_env_wksp")

### GET ESTIMATES:

preds.store<-list()

head(gds_upd,3); dim(gds_upd)

for (i in 1:nrow(best.mods)){

resp.thisrun<-as.character(best.mods$response[i])
pred.thisrun<-as.character(best.mods$predictor[i])

range.thisrun<-NULL

mod.thisrun<-NULL

type.thisrun<-as.character(best.mods$model[i])

if(type.thisrun=="add") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)+native,data=gds_upd)
if(type.thisrun=="int") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)*native,data=gds_upd)

summary(mod.thisrun)

nd.thisrun<-data.frame(native=c(rep("native",50),rep("non_native",50)),pred=seq(min(gds_upd[,pred.thisrun],na.rm=T),max(gds_upd[,pred.thisrun],na.rm=T),length.out=50))
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

nd2<-data.frame(native=c(rep("native",1),rep("non_native",1)),pred=mean(gds_upd[,pred.thisrun],na.rm=T))
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

best.mods # models to plot
preds.store # environmental estimates

plot.data<-best.mods
plot.data$resp_lab<-c(rep("Population growth rate",1),rep("log(Reprod. effort / m2)",1),rep("Neutral genetic diversity",2),rep("Adaptive genetic diversity",2))
plot.data$pred_lab<-c(rep("Mean temperature",2),rep("Temperature seasonality",1),rep("Mean precipitation (mm)",1),rep("Temperature seasonality",1),rep("Mean precipitation (mm)",1))
plot.data$space_after<-c(rep("no",6))
plot.data$lett<-letters[1:6]

quartz(title="",width=7.5,height=9, dpi=70, pointsize=20)
par(mfrow=c(3,2),mgp=c(2.7,0.9,0),oma=c(0,0,0.5,6.5),mar=c(4.5,4,1.5,1.7))

for (i in 1:nrow(plot.data)){

resp.thisrun<-as.character(plot.data$response[i])
pred.thisrun<-as.character(plot.data$predictor[i])
est.thisrun<-preds.store[[i]]

# Scaled data:
scdat.thisrun<-gds_upd

# Unscaled data:
if(gregexpr("_sc",pred.thisrun)[[1]][1]>0) unsc_x<-substr(pred.thisrun,1,nchar(pred.thisrun)-3) else unsc_x<-pred.thisrun

unscdat.thisrun<-gds_upd[,c(resp.thisrun,pred.thisrun)]
head(unscdat.thisrun,2)

raw.scaled<-scdat.thisrun[,c(which(colnames(scdat.thisrun) %in% c("site_code","native", pred.thisrun)),which(colnames(scdat.thisrun)==resp.thisrun))]
raw.unscaled<-scdat.thisrun[,c(which(colnames(scdat.thisrun) %in% c("site_code","native",unsc_x)),which(colnames(scdat.thisrun)==resp.thisrun))]
head(raw.scaled); dim(raw.scaled)
head(raw.unscaled); dim(raw.unscaled)

yl.now<-c(min(c(min(est.thisrun$lci),min(raw.scaled[,resp.thisrun],na.rm=T))),max(c(max(est.thisrun$uci),max(raw.scaled[,resp.thisrun],na.rm=T))))

head(est.thisrun)

plot(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],ylim=yl.now,type="n",bty="l",xlab="",ylab="",xaxt="n", las=1,cex.axis=0.8)

ylab.now<-plot.data$resp_lab[i]

if(resp.thisrun!="LOG_ros_reprod_m2") title(ylab=ylab.now) else title(ylab=bquote("log(reprod. effort/m"^2*")"))

if(plot.data$mod_name[i]!="LOG_ros_reprod_m2_sp_sc"){
pg.ci("predictor","est.thisrun",x.subset="native",colour=rgb(0,0,0,0.1))
lines(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],lty=1)
lines(est.thisrun[est.thisrun$native=="non_native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="non_native"],lty=2,col="red")
}

points(raw.scaled[raw.scaled$native=="native",pred.thisrun],raw.scaled[raw.scaled$native=="native",resp.thisrun],col="black",cex=0.5, pch=20)
points(raw.scaled[raw.scaled$native=="non_native",pred.thisrun],raw.scaled[raw.scaled$native=="non_native",resp.thisrun],col="red",cex=0.5, pch=20)

xax<-seq(round(min(raw.scaled[,pred.thisrun]),0),round(max(raw.scaled[,pred.thisrun]),0),length.out=5)
head(raw.unscaled)
xaxl<-round(seq(round(min(raw.unscaled[,unsc_x]),0),round(max(raw.unscaled[,unsc_x]),0),length.out=5),0)

if(unsc_x=="ap") xaxl<-round(xaxl,-2)

if(unsc_x=="mt") 
{
new.axis<-seq(0,round(max(raw.unscaled[,unsc_x]),0),length.out=5)

cantthink<-c(raw.unscaled[,unsc_x])
new.pos<-predict.lm(lm(raw.scaled[,pred.thisrun]~ cantthink),newdata=data.frame(cantthink=new.axis))

axis(side=1, at=new.pos, labels=new.axis, cex.axis=0.8)

} # close new axis

if(unsc_x!="mt") axis(side=1, at=xax, labels=xaxl, cex.axis=0.8)

xlab.thirun<-as.character(plot.data$pred_lab[i])
if(length(grep("emperatu",xlab.thirun))==0) title(xlab= xlab.thirun,mgp=c(2.2,1,0)) else title(xlab=bquote(.(xlab.thirun)*" ("*degree*"C)"),mgp=c(2.2,1,0))

mtext(paste("(",plot.data$lett[i],")",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.8,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/7))

if(plot.data$space_after[i]=="yes") blankplot()

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
supp.mods<-supp.mods[-which(supp.mods$model=="pred"),]
supp.mods<-rbind(supp.mods,sp_reprod)
supp.mods<-tidy.df(supp.mods)
supp.mods

# arrange for plotting:
supp.mods<-supp.mods[c(2,9,1,3,4,7,8,5,6),]
supp.mods<-tidy.df(supp.mods)
supp.mods

### GET ESTIMATES:

preds.store.supp<-list()

for (i in 1:nrow(supp.mods)){

resp.thisrun<-as.character(supp.mods$response[i])
pred.thisrun<-as.character(supp.mods$predictor[i])

range.thisrun<-NULL

mod.thisrun<-NULL

type.thisrun<-as.character(supp.mods$model[i])

if(type.thisrun=="add") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)+native,data=gds_upd)
if(type.thisrun=="int") mod.thisrun<-lm(get(resp.thisrun)~get(pred.thisrun)*native,data=gds_upd)

summary(mod.thisrun)

nd.thisrun<-data.frame(native=c(rep("native",50),rep("non_native",50)),pred=seq(min(gds_upd[,pred.thisrun],na.rm=T),max(gds_upd[,pred.thisrun],na.rm=T),length.out=50))
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

nd2<-data.frame(native=c(rep("native",1),rep("non_native",1)),pred=mean(gds_upd[,pred.thisrun],na.rm=T))
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

supp.mods # models to plot
preds.store.supp 

pd.supp<-supp.mods
pd.supp$resp_lab<-c(rep("log(reprod. effort/m2)",2),rep("Population growth rate",1),rep("Neutral genetic diversity",4),rep("Adaptive genetic diversity",2))
pd.supp$pred_lab<-c(rep("Mean temperature",1),rep("Precipitation seasonality",1),rep("Mean temperature",1),rep("Temperature seasonality",1),rep("Mean precipitation (mm)",1),rep("Density",1),rep("log(reprod. effort/m2)",1),rep("Mean temperature",1),rep("Mean precipitation (mm)",1))
pd.supp$space_after<-c(rep("no",2),"yes",rep("no",6))
pd.supp$lett<-letters[1:9]

quartz(title="",width=12.5,height=9.5, dpi=70, pointsize=21)
par(mfrow=c(3,4),mgp=c(2.7,0.9,0),oma=c(0,0,0.5,1),mar=c(4.5,4,1.5,1.7))

for (i in 1:nrow(pd.supp)){

resp.thisrun<-as.character(pd.supp$response[i])
pred.thisrun<-as.character(pd.supp$predictor[i])
est.thisrun<-preds.store.supp[[i]]
head(est.thisrun)

# Scaled data:
scdat.thisrun<-gds_upd

# Unscaled data:
if(gregexpr("_sc",pred.thisrun)[[1]][1]>0) unsc_x<-substr(pred.thisrun,1,nchar(pred.thisrun)-3) else unsc_x<-pred.thisrun

unscdat.thisrun<-gds_upd[,c(resp.thisrun,pred.thisrun)]
head(unscdat.thisrun,3)

raw.scaled<-scdat.thisrun[,c(which(colnames(scdat.thisrun) %in% c("site_code","native", pred.thisrun)),which(colnames(scdat.thisrun)==resp.thisrun))]
raw.unscaled<-scdat.thisrun[,c(which(colnames(scdat.thisrun) %in% c("site_code","native",unsc_x)),which(colnames(scdat.thisrun)==resp.thisrun))]
head(raw.scaled); dim(raw.scaled)
head(raw.unscaled); dim(raw.unscaled)

head(est.thisrun,3); dim(est.thisrun)

# For sp, only plot estimates within the range of the data:
if(pred.thisrun=="sp_sc"){
max.sp<-max(raw.scaled[which(raw.scaled$native=="native"),pred.thisrun])

est.thisrun<-est.thisrun[-which(est.thisrun$predictor[est.thisrun$native=="native"]>max.sp),]
est.thisrun<-tidy.df(est.thisrun)

}

if(length(grep("LOG_Y0_ros_m2",pred.thisrun))>0) raw.unscaled[,unsc_x]<-exp(raw.unscaled[,unsc_x])

yl.now<-c(min(c(min(est.thisrun$lci),min(raw.scaled[,resp.thisrun],na.rm=T))),max(c(max(est.thisrun$uci),max(raw.scaled[,resp.thisrun],na.rm=T))))

head(est.thisrun)

plot(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],ylim=yl.now,xlim=c(min(est.thisrun$predictor),max(est.thisrun$predictor)),type="n",bty="l",xlab="",ylab=pd.supp$resp_lab[i],xaxt="n", las=1,cex.axis=0.8)

pg.ci("predictor","est.thisrun",x.subset="native",colour=rgb(0,0,0,0.1))
lines(est.thisrun[est.thisrun$native=="native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="native"],lty=1)
lines(est.thisrun[est.thisrun$native=="non_native",which(colnames(est.thisrun)=="predictor")],est.thisrun$fit[est.thisrun$native=="non_native"],lty=2,col="red")

head(raw.scaled)
points(raw.scaled[raw.scaled$native=="native",pred.thisrun],raw.scaled[raw.scaled$native=="native",resp.thisrun],col="black",cex=0.5, pch=20)
points(raw.scaled[raw.scaled$native=="non_native",pred.thisrun],raw.scaled[raw.scaled$native=="non_native",resp.thisrun],col="red",cex=0.5, pch=20)

xax<-seq(round(min(raw.scaled[,pred.thisrun]),0),round(max(raw.scaled[,pred.thisrun]),0),length.out=5)
head(raw.unscaled)
xaxl<-round(seq(round(min(raw.unscaled[,unsc_x]),0),round(max(raw.unscaled[,unsc_x]),0),length.out=5),0)

if(unsc_x=="ap") xaxl<-round(xaxl,-2)
if(unsc_x=="LOG_Y0_ros_m2") xaxl<-round(xaxl,-1)
if(unsc_x=="sp") xaxl<-round(xaxl,-1)
if(unsc_x=="Y0_ros_m2") xaxl<-round(xaxl,-1)

if(unsc_x=="mt") 
{
new.axis<-seq(0,round(max(raw.unscaled[,unsc_x]),0),length.out=5)

ru_dat<-c(raw.unscaled[,unsc_x])
new.pos<-predict.lm(lm(raw.scaled[,pred.thisrun]~ ru_dat),newdata=data.frame(ru_dat=new.axis))

axis(side=1, at=new.pos, labels=new.axis, cex.axis=0.75)

} # close new axis

if(unsc_x!="mt") axis(side=1, at=xax, labels=xaxl, cex.axis=0.75)

xlab.thirun<-as.character(pd.supp$pred_lab[i])
if(length(grep("emperatu",xlab.thirun))==0) title(xlab=xlab.thirun,mgp=c(2.2,1,0)) else title(xlab=bquote(.(xlab.thirun)*" ("*degree*"C)"),mgp=c(2.2,1,0))

if(pd.supp$model[i]=="int") mtext(paste("(",pd.supp$lett[i],")    ","interaction",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.6,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/7)) else mtext(paste("(",pd.supp$lett[i],")    ","no interaction",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.6,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/7))

# Add general legend:
if(i==3) {
par(xpd=NA)
blankplot()

legend(0,10,legend=c("native","non native"),lty=c(1,2),bty="n",col=c("black","red"))
legend(1,10,legend=c("",""),pch=c(20),col=c("black","red"),bty="n")
par(xpd=F)
} # close if i==1

} # close plots i

# -------------------- #
# CORRELATION plots:
# -------------------- #

# mt, st, ap

env.testdf<-gds_upd[,c("mt","st","ap")]
head(env.testdf)

env.cortest<-cor(env.testdf,use="na.or.complete")
env.cornames<-combn(colnames(env.testdf),2)

axis.names<-data.frame(cornames=unique(paste(env.cornames)),newnames=c("Mean temperature","Temperature seasonality","Mean precipitation (mm)"))

quartz(title="",width=6.5,height=6.5, dpi=90, pointsize=16)
par(mfrow=c(2,2),mgp=c(2.7,0.9,0),oma=c(0,0,0.5,1),mar=c(4.5,4,1.5,1.7),xpd=F)

for (i in 1:length(env.testdf)){

v1<-env.cornames[1,i]
v2<-env.cornames[2,i]

col1<-which(colnames(gds_upd)==v1)
col2<-which(colnames(gds_upd)==v2)

v1<-as.character(axis.names$newnames[which(axis.names$cornames==v1)])
v2<-as.character(axis.names$newnames[which(axis.names$cornames==v2)])

plot(gds_upd[gds_upd$native=="native",col1], gds_upd[gds_upd$native=="native", col2],type="p",bty="o",xlab="",ylab="", las=1,cex.axis=0.9, col="black",cex=0.5, pch=20,xlim=c(min(gds_upd[,col1]),max(gds_upd[,col1])), ylim=c(min(gds_upd[, col2]),max(gds_upd[, col2])))
points(gds_upd[gds_upd$native=="non_native",col1], gds_upd[gds_upd$native=="non_native",col2],col="red",pch=20, cex=0.5)

if(length(grep("emperatu",v1))==0) title(xlab=v1,mgp=c(2.2,1,0)) else title(xlab=bquote(.(v1)*" ("*degree*"C)"),mgp=c(2.2,1,0))
if(length(grep("emperatu",v2))==0) title(ylab=v2,mgp=c(2.8,1,0)) else title(ylab=bquote(.(v2)*" ("*degree*"C)"),mgp=c(2.2,1,0))

mtext(paste("(",letters[i],")",sep=""),side=2,las=1,line=1.6,adj=0,cex=0.9,at=par("usr")[4]+((par("usr")[4]-par("usr")[3])/10))

st_mt<-cor.test(gds_upd[,col1], gds_upd[,col2])
text(par("usr")[2]-((par("usr")[2]/100)*4),(par("usr")[4]/6)*5.5,bquote(italic("r")*" = "*.(round(st_mt$estimate,2))),cex=0.7,adj=1)
text(par("usr")[2]-((par("usr")[2]/100)*4),(par("usr")[4]/6)*5,bquote(italic("p")*" = "*.(round(st_mt$p.value,2))),cex=0.7,adj=1)

# Add general legend:
par(xpd=NA)
if(i==1) legend(par("usr")[2]+((par("usr")[2]/100)*20),par("usr")[4]+((par("usr")[4]/10)*0.05),legend=c("native","non native"),pch=c(20),col=c("black","red"),bty="n")
par(xpd=F)

if(i==1) blankplot()

} # close i

} # close gen div



























