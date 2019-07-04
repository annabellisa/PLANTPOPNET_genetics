
# author: Annabel Smith

# V8 was updated for making custom colour palettes and organising colours so that the colours are consistent for the major groups across multiple K plots 

str_plot_V8<-function(K,cluster.data,site.data,las.opt,yaxs.loc,col.pal,site.lab,add.ccode=NULL,...){

# put cultivar and outgroup at the end:
cluster.data$native<-factor(cluster.data$native,levels=c("native","non_native","cultivar","outgroup"))
cluster.data<-cluster.data[order(cluster.data$native),]
cluster.data<-tidy.df(cluster.data)
head(cluster.data)

# Organise columns:
cluster.data$max_p<-apply(cluster.data[,3:(2+K)],1,function(x) which(x==max(x)))

cluster.data<-cbind(cluster.data[,1:2],cluster.data[,3:(K+2)][,unique(cluster.data$max_p)],cluster.data[,(K+3):length(cluster.data)])

# make assignment data matrix:
dat.thisrun<-apply(cluster.data[grep("assig",colnames(cluster.data))],1,rbind)
head(dat.thisrun)

# make site data:
sdat.thisrun<-cluster.data[,c("site_code","c_code",colnames(cluster.data)[which(colnames(cluster.data)=="native"):length(cluster.data)])]

# Set colour scheme:
if(col.pal %in% rownames(brewer.pal.info)) cols.thisrun<-brewer.pal(K,col.pal) else cols.thisrun<-get(col.pal)[1:K]

# Plot:
par(mar=c(5,5,1,0),mgp=c(1,yaxs.loc+1,yaxs.loc),xpd=F)
p1<-barplot(dat.thisrun,cex.axis=1.5,col=cols.thisrun,space=0,border=cols.thisrun,xaxt="n",las=las.opt,ylab="Probability \nof assignment",cex.lab=1.5)

# Add lines for individuals and sites:
arrows(which(!duplicated(sdat.thisrun$site_code))-1,-0.08,which(!duplicated(sdat.thisrun$site_code))-1,0.995,code=0,lwd=1)
arrows(nrow(sdat.thisrun),-0.08,nrow(sdat.thisrun),0.995,code=0,lwd=1)
arrows(1:nrow(sdat.thisrun),0,1:nrow(sdat.thisrun),0.995,code=0,lwd=0.1)

# Draw a box:
arrows(0,1,nrow(sdat.thisrun),1,code=0,lwd=1)
arrows(0,0,nrow(sdat.thisrun),0,code=0,lwd=1)

# Add labels for main sites
if(site.lab=="site_code") axis(1,which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))+3,labels=toupper(sdat.thisrun$site_code[which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))]),tick=F,line=2.5,las=las.opt)

if(site.lab=="c_code") axis(1,which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))+3,labels=toupper(sdat.thisrun$c_code[which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))]),tick=F,line=2.5,las=las.opt)

if(add.ccode==T){

axis(3,which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))+3,labels=toupper(sdat.thisrun$c_code[sdat.thisrun$native!="outgroup"][which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))]),tick=F,line=2.5,las=las.opt)

} # close add.ccode

# add lines for regions
par(xpd=NA)

arrows(head(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1),-0.25,tail(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1),-0.25,lwd=2,code=0)
mtext("Europe (north to south)",at=which(sdat.thisrun$region=="Europe")[length(which(sdat.thisrun$region=="Europe"))/2],side=1,line=4, cex=1.5)

arrows(head(c(which(sdat.thisrun$c_code=="AU"),which(sdat.thisrun$c_code=="NZ")),1),-0.25,tail(c(which(sdat.thisrun$c_code=="AU"),which(sdat.thisrun$c_code=="NZ")),1),-0.25,lwd=2,code=0)

aunz_lines<-c(which(sdat.thisrun$c_code=="AU"),which(sdat.thisrun$c_code=="NZ"))[order(c(which(sdat.thisrun$c_code=="AU"),which(sdat.thisrun$c_code=="NZ")))]
mtext("Australia / NZ",at=aunz_lines[length(aunz_lines)/2],side=1,line=4, cex=1.5)

head(sdat.thisrun)
axis(1,which(sdat.thisrun$site_code=="TK")[round(length(which(sdat.thisrun$site_code=="TK"))/2,0)],labels="Asia",tick=F,line=5,las=2,cex.axis=1.5)

arrows(head(which(sdat.thisrun$region=="Nth_America"),1),-0.35,tail(which(sdat.thisrun$region=="Nth_America"),1),-0.35,lwd=2,code=0)
mtext("North America",at=which(sdat.thisrun$region=="Nth_America")[length(which(sdat.thisrun$region=="Nth_America"))/2],side=1,line=5, cex=1.5)

axis(1,which(sdat.thisrun$site_code=="SA")[round(length(which(sdat.thisrun$site_code=="SA"))/2,0)],labels="Africa",tick=F,line=5,las=2,cex.axis=1.5)
axis(1,grep("OG",sdat.thisrun$site_code)[round(length(grep("OG",sdat.thisrun$site_code))/2,0)],labels="Outgroups",tick=F,line=2.5,las=2,cex.axis=1.5)

arrows(head(which(sdat.thisrun$native=="cultivar"),1),-0.25,tail(which(sdat.thisrun$native=="cultivar"),1),-0.25,lwd=2,code=0)
mtext("Cultivar",at=which(sdat.thisrun$native=="cultivar")[length(which(sdat.thisrun$native=="cultivar"))/2],side=1,line=4, cex=1.5)

# Add labels for native / non-native:

arrows(head(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1),1.05,tail(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1)-1,1.05,lwd=2,code=3,angle=90,length=0.05)
mtext("Native",at=which(sdat.thisrun$region=="Europe")[length(which(sdat.thisrun$region=="Europe"))/2],side=3,line=1, cex=1.5)

arrows(head(which(sdat.thisrun$region=="Nth_America"),1)+1,1.05,tail(which(sdat.thisrun$site_code=="SA"),1),1.05,lwd=2,code=3,angle=90,length=0.05)
mtext("Non-native",at=head(which(sdat.thisrun$region=="Nth_America"),1)+((tail(which(sdat.thisrun$site_code=="SA"),1)-head(which(sdat.thisrun$region=="Nth_America"),1))/2),side=3,line=1, cex=1.5)

par(xpd=F)

} # close structure plot function V8

# V7 was updated to put site codes on the xaxis and country codes on the top axis

str_plot_V7<-function(K,cluster.data,site.data,las.opt,yaxs.loc,col.pal,site.lab,add.ccode=NULL,...){

# put cultivar and outgroup at the end:
cluster.data$native<-factor(cluster.data$native,levels=c("native","non_native","cultivar","outgroup"))
cluster.data<-cluster.data[order(cluster.data$native),]
cluster.data<-tidy.df(cluster.data)

# make assignment data matrix:
dat.thisrun<-apply(cluster.data[grep("assig",colnames(cluster.data))],1,rbind)

# make site data:
sdat.thisrun<-cluster.data[,c("site_code","c_code",colnames(cluster.data)[which(colnames(cluster.data)=="native"):length(cluster.data)])]

# Set colour scheme:
cols.thisrun<-brewer.pal(K,col.pal)

# Plot:
par(mar=c(5,5,1,0),mgp=c(1,yaxs.loc+1,yaxs.loc),xpd=F)
p1<-barplot(dat.thisrun,cex.axis=1.8,col=cols.thisrun,space=0,border=cols.thisrun,xaxt="n",las=las.opt,ylab="Probability \nof assignment",cex.lab=1.8)

# Add lines for individuals and sites:
arrows(which(!duplicated(sdat.thisrun$site_code))-1,-0.08,which(!duplicated(sdat.thisrun$site_code))-1,0.995,code=0,lwd=1)
arrows(nrow(sdat.thisrun),-0.08,nrow(sdat.thisrun),0.995,code=0,lwd=1)
arrows(1:nrow(sdat.thisrun),0,1:nrow(sdat.thisrun),0.995,code=0,lwd=0.1)

# Draw a box:
arrows(0,1,nrow(sdat.thisrun),1,code=0,lwd=1)
arrows(0,0,nrow(sdat.thisrun),0,code=0,lwd=1)

# Add labels for main sites
if(site.lab=="site_code") axis(1,which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))+3,labels=toupper(sdat.thisrun$site_code[which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))]),tick=F,line=2.5,las=las.opt)

if(site.lab=="c_code") axis(1,which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))+3,labels=toupper(sdat.thisrun$c_code[which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))]),tick=F,line=2.5,las=las.opt)

if(add.ccode==T){

axis(3,which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))+3,labels=toupper(sdat.thisrun$c_code[sdat.thisrun$native!="outgroup"][which(!duplicated(sdat.thisrun$site_code[sdat.thisrun$native!="outgroup"]))]),tick=F,line=2.5,las=las.opt)

} # close add.ccode

# add lines for regions
par(xpd=NA)

arrows(head(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1),-0.25,tail(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1),-0.25,lwd=2,code=0)
mtext("Europe (north to south)",at=which(sdat.thisrun$region=="Europe")[length(which(sdat.thisrun$region=="Europe"))/2],side=1,line=4.5, cex=1.8)

arrows(head(which(sdat.thisrun$region=="Australasia" & sdat.thisrun$native=="non_native"),1),-0.25,tail(which(sdat.thisrun$region=="Australasia" & sdat.thisrun$native=="non_native"),1),-0.25,lwd=2,code=0)
mtext("Australasia",at=which(sdat.thisrun$region=="Australasia")[length(which(sdat.thisrun$region=="Australasia" & sdat.thisrun$native=="non_native"))/2],side=1,line=4.5, cex=1.8)

arrows(head(which(sdat.thisrun$region=="Nth_America"),1),-0.35,tail(which(sdat.thisrun$region=="Nth_America"),1),-0.35,lwd=2,code=0)
mtext("North America",at=which(sdat.thisrun$region=="Nth_America")[length(which(sdat.thisrun$region=="Nth_America"))/2],side=1,line=5.5, cex=1.8)

axis(1,which(sdat.thisrun$site_code=="SA")[round(length(which(sdat.thisrun$site_code=="SA"))/2,0)],labels="Africa",tick=F,line=5,las=2,cex.axis=1.8)
axis(1,grep("OG",sdat.thisrun$site_code)[round(length(grep("OG",sdat.thisrun$site_code))/2,0)],labels="Outgroups",tick=F,line=2.5,las=2,cex.axis=1.8)

arrows(head(which(sdat.thisrun$native=="cultivar"),1),-0.25,tail(which(sdat.thisrun$native=="cultivar"),1),-0.25,lwd=2,code=0)
mtext("Cultivar",at=which(sdat.thisrun$native=="cultivar")[length(which(sdat.thisrun$native=="cultivar"))/2],side=1,line=4.5, cex=1.5)

# Add labels for native / non-native:

arrows(head(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1),1.05,tail(which(sdat.thisrun$region=="Europe" & sdat.thisrun$native=="native"),1)-1,1.05,lwd=2,code=3,angle=90,length=0.05)
mtext("Native",at=which(sdat.thisrun$region=="Europe")[length(which(sdat.thisrun$region=="Europe"))/2],side=3,line=1, cex=1.8)

arrows(head(which(sdat.thisrun$region=="Nth_America"),1)+1,1.05,tail(which(sdat.thisrun$site_code=="SA"),1),1.05,lwd=2,code=3,angle=90,length=0.05)
mtext("Non-native",at=head(which(sdat.thisrun$region=="Nth_America"),1)+((tail(which(sdat.thisrun$site_code=="SA"),1)-head(which(sdat.thisrun$region=="Nth_America"),1))/2),side=3,line=1, cex=1.8)

par(xpd=F)



} # close structure plot function V7







