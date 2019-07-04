
# ------------------------------------ #
# ----------- SUPPLEMENT 05  --------- #
# ------------------------------------ #

### process demography
### Author: Annabel Smith

# load functions:
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

#########################################
####  	     	  DATA: 	   		 ####
#########################################

dat_dir<-"../01_data"
dir(dat_dir)

dem<-read.table(paste(dat_dir,"demog.txt",sep="/"),header=T)
dem$pid<-paste(dem$transect,dem$plot,dem$plant_id,sep="_")
dem$site_yr<-as.factor(paste(dem$site_code,dem$s_year,sep="_"))
head(dem)

#########################################
####  	     	  DENSITY: 	   		 ####
#########################################

dens<-dem[,c("site_code","c_year","s_year","transect","plot","plant_id","survival","no_rosettes","pid")]
head(dens)

# Most have two, some have one and some have three:
dens_sum<-as.data.frame.matrix(table(dens$site_code, dens$c_year))
dens_sum<-data.frame(site_code=rownames(dens_sum),dens_sum)
dens_sum<-tidy.df(dens_sum)
colnames(dens_sum)[2:length(dens_sum)]<-unique(dens$c_year)[order(unique(dens$c_year))]

dens_sum$Y0_plant_m2<-NA
dens_sum$Y1_plant_m2<-NA
dens_sum$Y2_plant_m2<-NA
dens_sum$Y0_ros_m2<-NA
dens_sum$Y1_ros_m2<-NA
dens_sum$Y2_ros_m2<-NA
head(dens_sum)
head(dens)

# CHECK: where survival=="no" the number of rosettes should only ever be NA or 1:
table(dem[which(dem$survival=="no"),]$no_rosettes)

# Re check the number of rosettes is == to the number of rows for each alive plant (see code below): TRUE

# Add DENSITY:

sy<-levels(dem$site_yr)

for (i in 1:length(sy)){

sy.thisrun<-sy[i]
dat.thisrun<-dem[dem$site_yr==sy.thisrun,]
dat.thisrun<-tidy.df(dat.thisrun)
head(dat.thisrun)

yr.thisrun<-levels(dat.thisrun$s_year)

if (yr.thisrun=="Y0"){

n_plants<-length(levels(dat.thisrun$plant_id))
n_plots<-length(levels(dat.thisrun$plot))
site.thisrun<-levels(dat.thisrun$site_code)

# I excluded na.rm=T on purpose because NAs flag data errors:
if(length(which(duplicated(dat.thisrun$plant_id)))>0) ros.dat<-dat.thisrun$no_rosettes[-which(duplicated(dat.thisrun$plant_id))] else ros.dat<-dat.thisrun$no_rosettes
if(length(ros.dat)!=n_plants) stop("something is wrong")
n_ros<-sum(ros.dat)

# Plants per m2:
dens_sum[which(dens_sum==site.thisrun),grep(paste(yr.thisrun,"plant",sep="_"),colnames(dens_sum))]<-(n_plants/n_plots)*4

# Ros per m2:
dens_sum[which(dens_sum==site.thisrun),grep(paste(yr.thisrun,"ros",sep="_"),colnames(dens_sum))]<-(n_ros/n_plots)*4

} # close if Y0

if (yr.thisrun!="Y0"){

# Remove plants which didn't survive; everything else should be a survivor or a new plant:
if(length(which(dat.thisrun$survival=="no"))>0) surv_plant<-dat.thisrun[-which(dat.thisrun$survival=="no"),] else surv_plant<-dat.thisrun
surv_plant<-tidy.df(surv_plant)

n_plants<-length(levels(surv_plant$plant_id))
n_plots<-length(levels(dat.thisrun$plot))
site.thisrun<-levels(dat.thisrun$site_code)

# Here it's necessary to remove NAs:
if(length(which(duplicated(surv_plant$plant_id)))>0) ros.dat<-surv_plant$no_rosettes[-which(duplicated(surv_plant$plant_id))] else ros.dat<-surv_plant$no_rosettes
if(length(which(is.na(ros.dat)))>0) ros.dat<-ros.dat[-which(is.na(ros.dat))] else ros.dat<-ros.dat

# This line will detect data errors:
if(length(ros.dat)!=n_plants) stop("something is wrong")
n_ros<-sum(ros.dat)

# Plants per m2:
dens_sum[which(dens_sum==site.thisrun),grep(paste(yr.thisrun,"plant",sep="_"),colnames(dens_sum))]<-(n_plants/n_plots)*4

# Ros per m2:
dens_sum[which(dens_sum==site.thisrun),grep(paste(yr.thisrun,"ros",sep="_"),colnames(dens_sum))]<-(n_ros/n_plots)*4

} # close if not Y0

} # close for i

# Check that the right number of cols were added:
dens_sum[which(unlist(apply(dens_sum,1,function(x) length(which(x[2:5]>0))==length(which(!is.na(x[6:8])))))==F),]

# Plot density of plants against density of rosettes:

quartz("",12,4,pointsize=20,dpi=80)
par(mfrow=c(1,3), mar=c(5,5,1,1))

yrs<-c("Y0","Y1","Y2")

for (i in 1:length(yrs)){

yr.thisrun<-yrs[i]
x.thisrun<-grep(paste(yr.thisrun,"_plant_m2",sep=""),colnames(dens_sum))
y.thisrun<-grep(paste(yr.thisrun,"_ros_m2",sep=""),colnames(dens_sum))

plot(dens_sum[,x.thisrun], dens_sum[,y.thisrun], pch=20, las=1, xlab="Plants / m2", ylab="Rosettes / m2", xlim=c(0,400))

text(0,par("usr")[4]-(par("usr")[4]/6),paste(yr.thisrun,"\nR2 = ",round(cor(dens_sum[,x.thisrun], dens_sum[,y.thisrun], use="complete.obs"),2),sep=""),adj=0, col="cornflowerblue",font=2)

} # close for

# Calculate change in density (empirical population growth rate):

dens_sum$emp_pgr_plt<-NA
dens_sum$emp_pgr_ros<-NA
head(dens_sum)

sites<-as.character(dens_sum$site_code)

for (i in 1:length(sites)){

site.thisrun<-as.character(dens_sum$site_code[i])
data.thisrun<-dens_sum[dens_sum$site_code==site.thisrun,]

plt.thisrun<-data.thisrun[,grep("plant_m2",colnames(data.thisrun))]
if(length(which(is.na(plt.thisrun)))>0) plt.thisrun<-plt.thisrun[-which(is.na(plt.thisrun))]
if(length(plt.thisrun)>1) plt.thisrun<-plt.thisrun[1:2] else next

ros.thisrun<-data.thisrun[,grep("ros_m2",colnames(data.thisrun))]
if(length(which(is.na(ros.thisrun)))>0) ros.thisrun<-ros.thisrun[-which(is.na(ros.thisrun))]
if(length(ros.thisrun)>1) ros.thisrun<-ros.thisrun[1:2] else next

dens_sum$emp_pgr_plt[which(dens_sum$site_code==site.thisrun)]<-log(plt.thisrun[2]/plt.thisrun[1])
dens_sum$emp_pgr_ros[which(dens_sum$site_code==site.thisrun)]<-log(ros.thisrun[2]/ros.thisrun[1])

} # close for

#########################################
####  	       FECUNDITY: 	   		 ####
#########################################

# Add flowering col:
dem$flowering<-ifelse(dem$no_fl_stems>0,1,0)
head(dem,4)

# Add reproductive effort:
dem$reprod<-dem$no_fl_stems * dem$inflor_length

# Where no_fl_stems==0 make the reprod col 0
dem$reprod[dem$no_fl_stems==0]<-0
check.rows(dem[,c("site_code","s_year","pid","survival","no_rosettes","no_fl_stems","fl_stem_height","inflor_length","site_yr","flowering","reprod")])

# Size corrected reproductive effort:
head(dem,3)
dem$repr_sc<-dem$reprod/dem$leaf_length
# plot(dem$reprod, dem$repr_sc)

# remove outliers:
dem<-dem[-which(dem$repr_sc>40),]
dem<-tidy.df(dem)

head(dem)

# PLANT LEVEL (not rosette) FECUNDITY
# Y0 ONLY
head(dem,3)

# Rosette level data:
rsf<-dem[which(dem$s_year=="Y0"),]
rsf<-tidy.df(rsf)
rsf$site_pid<-paste(rsf$site_code, rsf$pid, sep="_")
head(rsf,3); dim(rsf)

# Plant level data:

# flowering = whole plant is flowering (1) or not (0):
plf<-aggregate(flowering~site_pid,FUN=sum,na.rm=T,data=rsf)
plf$flowering<-ifelse(plf$flowering==0,0,1)

# reprod = sum of total reprod effort across all rosettes in plant:
plrep<-aggregate(reprod~site_pid,FUN=sum,na.rm=T,data=rsf)
# reprod = sum of total size corrected reprod effort across all rosettes in plant:
plrepsc<-aggregate(repr_sc~site_pid,FUN=sum,na.rm=T,data=rsf)
head(plrepsc); dim(plrepsc)

head(rsf,3)
plfecund<-rsf[,c("site_code","c_year","s_year","transect","plot","plant_id","x_coord","y_coord","survival","no_rosettes","pid","site_yr","site_pid")]
plfecund<-plfecund[-which(duplicated(plfecund$site_pid)),]
plfecund<-tidy.df(plfecund)
head(plfecund); dim(plfecund)

# Check all names occur in each of the three new data frames:
table(plrepsc$site_pid %in% plfecund$site_pid)

# Then merge:
plfecund<-merge(plfecund,plf,by="site_pid",all.x=T, all.y=F)
plfecund<-merge(plfecund,plrep,by="site_pid",all.x=T, all.y=F)
plfecund<-merge(plfecund,plrepsc,by="site_pid",all.x=T, all.y=F)

head(plfecund); dim(plfecund)
range(plfecund$flowering)
# there are 78 and 86 NAs in the reprod and repr_sc cols respectively, which look like cases of no record for the flowering stem, perhaps because they were broken etc. The flowering col. can still be used in the prop flowering calc, and the NAs can be removed in the means for their others:
plfecund[(which(is.na(plfecund$reprod))[1:10]),]

# Calculate site x yr fecundity:
# Done only for Y0 at this stage

# Remove extreme outliers (e.g. over 2000)? Haven't done this as they look feasible
# 1. proportion flowering
# 2. average reproductive effort = mean(reprod)
# 2. density of reproductive effort = reprod / number of plots 

dens_sum$rs_prop_fl<-NA
dens_sum$rs_mean_reprod<-NA
dens_sum$rs_mean_reprod_sc<-NA

dens_sum$rs_reprod_m2<-NA
dens_sum$rs_reprod_sc_m2<-NA

dens_sum$pl_prop_fl<-NA
dens_sum$pl_mean_reprod<-NA
dens_sum$pl_mean_reprod_sc<-NA

dens_sum$pl_reprod_m2<-NA
dens_sum$pl_reprod_sc_m2<-NA

sites<-levels(dens_sum$site_code)

# Rosette level data (ALL YEARS, reduced in loop):
head(dem,3); dim(dem)

# Plant level data (Y0 only):
head(plfecund,3); dim(plfecund)

for (i in 1:length(sites)){

site.thisrun<-sites[i]

# Rosette level data:
dat.thisrun<-dem[dem$site_code== site.thisrun,]
dat.thisrun<-dat.thisrun[dat.thisrun$s_year=="Y0",]
dat.thisrun<-tidy.df(dat.thisrun)
head(dat.thisrun)

if(levels(dat.thisrun$s_year)!="Y0") stop ("something is wrong")

# Plant level data:
pl.thisrun<-plfecund[plfecund$site_code==site.thisrun,]
pl.thisrun<-tidy.df(pl.thisrun)
head(pl.thisrun)

n_plots<-length(levels(dat.thisrun$plot))

# this is wrong for reproductive density: 
# dens_sum$reprod_m2_old[i]<-(mean(dat.thisrun$reprod, na.rm=T)/n_plots)*4

# -- ROSETTE LEVEL:
dens_sum$rs_prop_fl[i]<-mean(dat.thisrun$flowering,na.rm=T)
dens_sum$rs_mean_reprod[i]<-mean(dat.thisrun$reprod, na.rm=T)
dens_sum$rs_mean_reprod_sc[i]<-mean(dat.thisrun$repr_sc, na.rm=T)

dens_sum$rs_reprod_m2[i]<-(sum(dat.thisrun$reprod, na.rm=T)/n_plots)*4
dens_sum$rs_reprod_sc_m2[i]<-(sum(dat.thisrun$repr_sc, na.rm=T)/n_plots)*4

# -- PLANT LEVEL:
head(pl.thisrun)

dens_sum$pl_prop_fl[i]<-mean(pl.thisrun$flowering,na.rm=T)
dens_sum$pl_mean_reprod[i]<-mean(pl.thisrun$reprod, na.rm=T)
dens_sum$pl_mean_reprod_sc[i]<-mean(pl.thisrun$repr_sc, na.rm=T)

dens_sum$pl_reprod_m2[i]<-(sum(pl.thisrun$reprod, na.rm=T)/n_plots)*4
dens_sum$pl_reprod_sc_m2[i]<-(sum(pl.thisrun$repr_sc, na.rm=T)/n_plots)*4

} # close i

head(dens_sum,3); dim(dens_sum)

# Merge with site data:
sdat<-read.table(paste(dat_dir,"site_data.txt",sep="/"),header=T)
sdat2<-sdat[,1:13]
head(sdat2,2)

# Update names:
dens_sum$site_code<-as.character(dens_sum$site_code)
dens_sum$site_code[dens_sum$site_code=="LK"]<-"LK1"
dens_sum$site_code<-as.factor(dens_sum$site_code)
dens_sum$emp_pgr_plt<-unlist(dens_sum$emp_pgr_plt)
dens_sum$emp_pgr_ros<-unlist(dens_sum$emp_pgr_ros)
dens_sum<-tidy.df(dens_sum)

# Check all names are present in site_data (this should be zero):
length(which(dens_sum$site_code %in% sdat2$site_code==F))

ds<-dens_sum[,c(1, which(colnames(dens_sum)=="Y0_plant_m2"):length(dens_sum))]

head(sdat)
head(ds)

sdt<-merge(sdat2, ds, by="site_code", all.x=T, all.y=F)
head(sdt,3)

# write.table(sdt, "sdt.txt", sep="\t", row.names=F, quote=F)








































