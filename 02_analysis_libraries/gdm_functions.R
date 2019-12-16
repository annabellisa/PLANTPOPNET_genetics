
# Author: Annabel Smith

# Format data for GDM and plot outputs. Based on functions in gdm package but modified to provide more control.

format_gdm<-function(orig_data,vars.toinclude,enviro.index){

# orig_data is a data frame including pairwise population statistics (fst, lat and lon for each site in the pair, and enviro variables for each site in each pair)
# vars.toinclude is a character vector of the variables in orig_data (colnames(orig_data)) that you want to include in the analysis. Enter these WITHOUT the 1 and 2 which designate sites within pairs. 
# enviro.index specifies which column the environmental variables start on, excluding fsts and lats and lons. 

gdm_dat<-data.frame(distance=orig_data$fst,weights=1,orig_data[,c(which(colnames(orig_data) %in% c("lat1","lon1","lat2","lon2")),which(colnames(orig_data) %in% colnames(orig_data)[enviro.index:length(colnames(orig_data))][unlist(sapply(vars.toinclude,function(x)grep(x,colnames(orig_data)[enviro.index:length(colnames(orig_data))])))]))])

# make native numeric:
if (length(grep("nat1",colnames(gdm_dat)))>0) gdm_dat$nat1<-as.numeric(ifelse(gdm_dat$nat1=="native",0,1))
if (length(grep("nat2",colnames(gdm_dat)))>0) gdm_dat$nat2<-as.numeric(ifelse(gdm_dat$nat2=="native",0,1))

# Re-arrange:
gdm_dat<-gdm_dat[,c(1:6,which(colnames(gdm_dat) %in% colnames(gdm_dat)[7:length(gdm_dat)][grep("1",colnames(gdm_dat)[7:length(gdm_dat)])]),which(colnames(gdm_dat) %in% colnames(gdm_dat)[7:length(gdm_dat)][grep("2",colnames(gdm_dat)[7:length(gdm_dat)])]))]

# Prefix all vars with s1. and s2.
colnames(gdm_dat)[grep("1",colnames(gdm_dat))]<-paste("s1.",colnames(gdm_dat)[grep("1",colnames(gdm_dat))],sep="")
colnames(gdm_dat)[grep("2",colnames(gdm_dat))]<-paste("s2.",colnames(gdm_dat)[grep("2",colnames(gdm_dat))],sep="")

# Remove trailing number:
colnames(gdm_dat)[3:length(colnames(gdm_dat))]<-substr(colnames(gdm_dat)[3:length(colnames(gdm_dat))],1,nchar(colnames(gdm_dat)[3:length(colnames(gdm_dat))])-1)

# Change names of coordinate cols:
colnames(gdm_dat)[3:6]<-c("s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord")
head(gdm_dat,3)

return(gdm_dat)

} # close format gdm

# simplified version for SI:
plot_gdm_simp<-function(spline_dat,gdm_fit,gdm_xlabs,gdm_test,adjust_p=NULL){

if(adjust_p==T) gdm_ps<-p.adjust(gdm_test[[3]][,1],method="bonferroni") else gdm_ps<-gdm_test[[3]][,1]

for (i in 1:(length(spline_dat$x[1,]))){

spline.thisrun<-colnames(spline_dat$x)[i]
if(gdm_ps[i]==0) p.prefix<-"P < " else p.prefix<-"P = "

if(spline.thisrun!="nat"){
x.thisrun<-spline_dat$x[,i]
y.thisrun<-spline_dat$y[,i]
plot(x.thisrun,y.thisrun,ylim=c(range(spline_dat$y)[1],range(spline_dat$y)[2]),cex=0.5,pch=20,xlab="",ylab=expression("Partial genetic distance ("*italic("F")[ST]*")"),type="l")

xlab.thisrun<-as.character(gdm_xlabs$newname[i])

if(length(grep("emperature",xlab.thisrun))==0) title(xlab=xlab.thisrun) else title(xlab=bquote(.(xlab.thisrun)*" ("*degree*"C)"))

text(min(x.thisrun),max(spline_dat$y)-0.005,paste("Max. spline = ",round(max(y.thisrun),3),sep=""),adj=0,col="red")

if(max(spline_dat$y)<0.1) ylim.ppos<-(max(spline_dat$y)/3.8)*3 else ylim.ppos<-(max(spline_dat$y)/3.5)*3 

if (gdm_ps[i]==0 & gdm_ps[i]!=-9999) text(min(x.thisrun),ylim.ppos,paste(p.prefix,"0.001",sep=""),adj=0,col="red") 

if (gdm_ps[i]!=0  & gdm_ps[i]!=-9999) text(min(x.thisrun),ylim.ppos,paste(p.prefix,round(gdm_ps[i],3),sep=""),adj=0,col="red")

if (gdm_ps[i]==-9999) text(min(x.thisrun),ylim.ppos,"P > 0.9",adj=0,col="red")

} # close spline not native

if(spline.thisrun=="nat"){
x.thisrun<-spline_dat$x[,i]
y.thisrun<-spline_dat$y[,i]

plot(c(0,1),c(head(y.thisrun,1),tail(y.thisrun,1)),lwd=2,xlab=gdm_xlabs$newname[i], ylab=expression("Partial genetic distance ("*italic("F")[ST]*")"),xlim=c(-0.5,1.5),ylim=c(0,max(spline_dat$y)),pch=20,xaxt="n")
axis(side=1,at=c(0,1),labels=c("same","different"))

text(-0.5,max(spline_dat$y)-0.005,paste("Max. spline = ",round(max(y.thisrun),3),sep=""),adj=0,col="red")

if(max(spline_dat$y)<0.1) ylim.ppos<-(max(spline_dat$y)/3.8)*3 else ylim.ppos<-(max(spline_dat$y)/3.5)*3 

if (gdm_ps[i]==0) text(-0.5,ylim.ppos,paste(p.prefix,"0.001",sep=""),adj=0,col="red") else text(-0.5,ylim.ppos,paste(p.prefix,gdm_ps[i],sep=""),adj=0,col="red")

} # close spline native

} # close plot all

} # close function

format_gdm_xlabs<-function(gdm_splines) return(gdm_xlabs_all[gdm_xlabs_all$spline %in% colnames(gdm_splines$x),])
