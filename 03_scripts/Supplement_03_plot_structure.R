
# ------------------------------------ #
# ----------- SUPPLEMENT 04  --------- #
# ------------------------------------ #

# Plot structure results:
### Author: Annabel Smith

# load functions:
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

library("RColorBrewer")
library(rworldmap)
library("plotrix")

#########################################
####  	       DATA SET UP:    		 ####
#########################################
{
# SET WORKING DIRECTORIES, FILES AND TIDY DATA:

# The project dir is the location of the structure input files:
proj_dir<-"/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSIS_RESULTS/STRUCTURE/STRUCTURE_DIR/DPlan18_filt6"
dir(proj_dir)

# The results dir contains the structure results, usually within proj_dir:
res_dir<-paste(proj_dir,"DPlan18_filt6_results",sep="/")
dir(res_dir)

file_name<-"DPlan18_filt6"

# this is usually .fam but can also be .str:
str_file<-paste(proj_dir,dir(proj_dir)[grep(".fam",dir(proj_dir))],sep="/")

# str_sites is site data that corresponds to genetic assignment probabilities:
str_sites<-read.table(str_file,colClasses=c(rep("character",2),rep("NULL",dim(read.table(str_file))[2]-2)),header=F)

dat_dir<-"../01_data"

# read main site data file:
sdat<-read.table(paste(dat_dir,"site_data.txt",sep="/"),header=T)
sdat$region<-factor(sdat$region,levels=c("Europe","Nth_America","Australasia"))
sdat<-sdat[sdat$genetics=="Y",]
sdat<-sdat[-which(sdat$site_code %in% c("RCH")),]
sdat<-sdat[-which(sdat$n_gt<7 & sdat$native!="outgroup"),]
sdat<-sdat[,c("site_code","n_gt","native","region","country","c_code","location","latitude","longitude")]
sdat<-tidy.df(sdat)
sdat<-sdat[order(sdat$native,sdat$region,sdat$latitude),]
sdat<-tidy.df(sdat)
head(sdat,2); dim(sdat)

} # close data set-up

#########################################
#  SET K AND ASSIGNMENT PROBABILITIES:  #
#########################################

# Set K and run this whole section

# SET K:
K<-6

{

# Get assignment probabilities:
# use out_dir or res_dir, depending on where results are:
assig<-read.table(paste(res_dir,dir(res_dir)[grep(paste("\\.",K,".meanQ",sep=""),dir(res_dir))],sep="/"),header=F)

# Combine site_data and assigment probs:
site_assig<-cbind(str_sites,assig)
site_assig$V1<-toupper(site_assig$V1)

# this is the order for .fam files
colnames(site_assig)<-c("site_code","indiv",paste("assig",1:K,sep=""))

sdat2<-sdat[,c("site_code","c_code","native","region","country","latitude","longitude")]
head(sdat2)
head(site_assig)

# should be zero:
length(unique(site_assig$site_code)[which(unique(site_assig$site_code) %in% as.character(unique(sdat2$site_code))!=T)])

all_dat<-merge(site_assig,sdat2,by="site_code",all.x=T,all.y=F)
all_dat<-all_dat[order(all_dat$native,all_dat$region,-all_dat$latitude,all_dat$indiv),]
all_dat<-tidy.df(all_dat)
head(all_dat)

} # close set K

#########################################
####  	       BAR PLOTS:    		 ####
#########################################

# Single barplot (for main document):

# The sequential palettes names are:
#  Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd

# The diverging palettes are 
# BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

switch.col<-brewer.pal(K,"Accent")[c(4,5,1,6,2,3)]

# Test colours:
quartz(title="Fig",width=4,height=4,dpi=80,pointsize=10)
barplot(matrix(data=6,nrow=6,ncol=1),col=brewer.pal(K,"Accent"),yaxt="n")
text(1,seq(1,35,length.out=6),c("non-native","cultivar","outgroup","nth central eur","atlantic","greece"))
text(0.5,seq(1,35,length.out=6),switch.col)

### MAIN PLOT:
quartz(title="Fig",width=16,height=4,dpi=80,pointsize=10)
par(mfrow=c(1,1),oma=c(3,0,2,0))
str_plot_V8(K,all_dat,sdat2,las.opt=2,yaxs.loc=-3,cex.axis=0.7,col.pal="switch.col",site.lab="c_code",add.ccode=F)

# Multiple barplots (for SI):
{

load("../04_workspaces/Supp_plot_structure")

# Data frames saved in workspace:
k6dat<-all_dat
k7dat<-all_dat
k8dat<-all_dat
k9dat<-all_dat
k10dat<-all_dat
k11dat<-all_dat
k12dat<-all_dat
k13dat<-all_dat

save.image("../04_workspaces/Supp_plot_structure")

head(all_dat)
head(k6dat)
head(k7dat)
head(k8dat)
head(k9dat)
head(k10dat)
head(k11dat)
head(k12dat)
head(k13dat)

# COLOURS can be entered as a brewer.pal name in quotes OR as a named vector of colours in quotes. 

# check colours:
# barplot(matrix(data=11,nrow=11,ncol=1),col=new.pal,yaxt="n")

new.pal<-c(rev(brewer.pal(8,"Accent")),brewer.pal(5,"BrBG"))

dat.toplot<-data.frame(K=6:13,data=paste("k",6:13,"dat",sep=""))

site.data<-sdat2
las.opt<-2
yaxs.loc<--3
site.lab<-"site_code"
add.ccode=T
col.pal<-"new.pal"

i<-8
cluster.thisrun<-get(as.character(dat.toplot$data[i]))
head(cluster.thisrun)
k.thisrun<-dat.toplot$K[i]
k.thisrun

quartz(title="Fig",width=16,height=4,dpi=80,pointsize=10)
par(mfrow=c(1,1),oma=c(4,0,1,0))
str_plot_V8(k.thisrun,cluster.thisrun,site.data,las.opt,yaxs.loc,col.pal,site.lab,add.ccode)

} # close SI

#########################################
##    		 MAP FOR PAPER:  	       ##
#########################################
{

# Remove outgroups and cultivars:
sdpie<-sdat[-c(which(sdat$native=="cultivar"),which(sdat$native=="outgroup")),]
sdpie<-tidy.df(sdpie)
head(sdpie,2); dim(sdpie)

# Load GBIF data:
pdat<-read.csv(paste(dat_dir,"GBIF_BIEN_Plantago.csv",sep="/"),header=T)
head(pdat)

# Make pie specific data
cluster.data<-all_dat

# Remove outgroups and cultivars:
cluster.data<-cluster.data[-c(which(cluster.data$native=="cultivar"),which(cluster.data$native=="outgroup")),]
cluster.data<-tidy.df(cluster.data)

# Organise columns:
cluster.data$max_p<-apply(cluster.data[,grep("assig",colnames(cluster.data))],1,function(x) which(x==max(x)))

cluster.data<-cbind(cluster.data[,1:2],cluster.data[,grep("assig",colnames(cluster.data))][,unique(cluster.data$max_p)],cluster.data[,(K+3):length(cluster.data)])
head(cluster.data)

# aggregate assigment probabilities per site:
agg_dat<-do.call(cbind,apply(cluster.data[,grep("assig",colnames(cluster.data))],2,function(x)aggregate(x~cluster.data$site_code,cluster.data[,2:length(cluster.data)],sum)))
agg_dat<-data.frame(site_code=agg_dat[,1],agg_dat[,colnames(agg_dat)[grep(".x",colnames(agg_dat))]])
head(agg_dat)

# check:
site_now<-sample(agg_dat$site_code,1)
colSums(cluster.data[which(cluster.data$site_code==site_now),grep("assig",colnames(cluster.data))])
agg_dat[agg_dat$site_code==site_now,]

# Add lats and longs to pie data:
latlong<-sdpie[,c("site_code","latitude","longitude")]
r1<-merge(agg_dat,latlong, by.x="site_code", by.y="site_code", all.x=T, all.y=F)
r1<-tidy.df(r1)
head(r1); dim(r1)

quartz("",11,6,dpi=70)
par(mar=c(2,2,2,2))
par(mar=c(0,0,0,0))
plot(countriesCoarseLessIslands,lwd=0.5,bg="white",col="grey80",border=F)
points(pdat$lon, pdat$lat, pch=20, col="grey50", cex=0.5)

# Plot pies:
for (i in 1:nrow(r1)){

data.thisrun<-r1[i,]

floating.pie(data.thisrun$longitude, data.thisrun$latitude, x=as.numeric(data.thisrun[grep("assig",colnames(data.thisrun))]),radius=3,  border=NA,col=switch.col)

} # close for i

# Plot EUR only:

quartz("",10,8,dpi=70)
par(mar=c(0,0,0,0))
plot(countriesCoarseLessIslands,lwd=0.5,bg="white",col="grey",border=F,xlim = c(-10,40),ylim = c(45, 45),asp = 1)
# need to add pdat below:
points(pdat$lon, pdat$lat, pch=20, col="grey50", cex=0.5)

for (i in 1:nrow(r1)){

data.thisrun<-r1[i,]

floating.pie(data.thisrun$longitude, data.thisrun$latitude, x=as.numeric(data.thisrun[grep("assig",colnames(data.thisrun))]),radius=2,  border=NA,col=switch.col)

} # close for i

} # close map for paper


































