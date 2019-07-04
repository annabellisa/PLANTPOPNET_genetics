
###########################
#  	SIMULATION FUNCTIONS  #
###########################

# Author: Annabel Smith

# Based on tutorials by Marco Andrello: https://sites.google.com/site/marcoandrello/metapopgen

std.err<-function(x) sd(x)/sqrt(length(x))

# Initial population structure. Currently uses default notation for m, n, z and nind, so no arguments needed:
pop.wrapper1<-function(){

# Genotypes number are assigned at random, by drawing the number of individuals of each genotype from a uniform distribution between [0,nind/m]
	
N1 <- array(NA,dim=c(m,n,z))	# Initial number of individuals

for (deme in 1:n){

for (age in 1:z){

for (genotype in 1:m){
N1[genotype,deme,age] <- nind
} # close for genotype

} # close for age

} # close for deme

return(N1)

} # close pop.wrapper1

# Functionise simulation, close to the island example code, with fst calcs:
sim.wrapper1<-function(nrep=1){

# Simulation
nrepl <- nrep

# Create the array for fst
fst <- array(NA,dim=c(nrepl,tmax))

# Run the simulations
for (i in 1:nrepl){

cat("Replicate",i,"\n") # Just to print on screen the advancement of the simulation

N<-sim.metapopgen.monoecious(input.type="array",N1 = N1, sigma = sigma,phi_M = phi_M, phi_F = phi_F,mu = mu, delta = delta, kappa0 = kappa0,T_max = tmax,save.res = F)
return(N)

flush.console()

# Calculate Fst
print("Calculating fst...")

for (t in 1:tmax){

fst[i,t] <- fst.global.monoecious(N,t)

} # close for calc fst

} # close simulation
} # close sim.wrapper1

# This only gives genotype frequencies, it doesn't do the fst calc:
sim.wrapper2<-function(nrep=1,...){

N.out<-list()

# Run the simulations
for (i in 1:nrep){

# cat("Replicate",i,"\n") # Just to print on screen the advancement of the simulation

N.out[[i]]<-sim.metapopgen.monoecious(input.type="array",N1 = N1, sigma = sigma,phi_M = phi_M, phi_F = phi_F,mu = mu, delta = delta, kappa0 = kappa0,T_max = tmax,save.res = F)

} # close i simulation reps

return(N.out)

} # close sim.wrapper2

# Calculate multi-locus heterozygosities. Currently works on the N object output from sim.wrapper2:
he.multi<-function(nrep=3,ndeme=n){

if (n>2) stop ("function not defined for > 2 demes")

he.out<-list()

for (i in 1:nrep){

N.thisrun<-N[[i]]

# deme 1:
d1A<-apply(N.thisrun[,1,1,],2,freq.all)[1,]
d1B<-apply(N.thisrun[,1,1,],2,freq.all)[2,]
d1He<-1-(d1A^2+d1B^2)

# deme 2:
d2A<-apply(N.thisrun[,2,1,],2,freq.all)[1,]
d2B<-apply(N.thisrun[,2,1,],2,freq.all)[2,]
d2He<-1-(d2A^2+d2B^2)

he.out[[i]]<-data.frame(He=c(d1He,d2He))

} # close i nrep

he.result<-data.frame(deme=c(rep("d1",length(d1He)),rep("d2",length(d2He))),do.call(cbind,he.out))

colnames(he.result)[2:length(he.result)]<-paste("rep",1:nrep,sep="")

return(he.result)

} # close he.multi

# Summarise deme-level population size. Currently works on the N object output from sim.wrapper2:
popsize.multi<-function(nrep=3,ndeme=n){

if (n>2) stop ("function not defined for > 2 demes")

psm.out<-list()

for (i in 1:nrep){

N.thisrun<-N[[i]]

d1N<-colSums(N.thisrun[,1,1,])
d2N<-colSums(N.thisrun[,2,1,])

psm.out[[i]]<-data.frame(N=c(d1N,d2N))

} # close i nrep

N.result<-data.frame(deme=c(rep("d1",length(d1N)),rep("d2",length(d2N))),do.call(cbind,psm.out))

colnames(N.result)[2:length(N.result)]<-paste("rep",1:nrep,sep="")

return(N.result)

} # close popsize.multi

# plot sim results for single locus:
n_he_plot<-function(d1sigma=NULL,d2sigma=NULL,Fphi=NULL,Mphi=NULL){
par(mar=c(4,4,3,4),mgp=c(2.5,1,0),new=F)
plot(1:100,N[1,1,1,],ylim=c(0,sum(N1)/2),type="n",col="gray80",xlab="Years",ylab="",bty="u",las=1)
title(main=bquote(sigma["d1"]~" = "~.(d1sigma)*"; "*sigma["d2"]~" = "~.(d2sigma)),line=2)
title(main=bquote(phi["F"]~" = "~.(Fphi)*"; "*phi["M"]~" = "~.(Mphi)),line=0)
title(ylab="Population size",mgp=c(3,1,0))
lines(1:100,colSums(N[,1,1,]))
lines(1:100,colSums(N[,2,1,]),col="red")
par(new=T)
plot(1:100,d1He,axes=F,type="l",xlab=NA,ylab="",ylim=c(0,1),las=1,lty=2)
axis(side=4,las=1)
mtext("He",side=4,ylab="xx",cex = par("cex.lab"),line=2.5,las=0)
lines(1:100,d2He,col="red",lty=2)
legend("bottomleft",legend=c("d1 N","d2 N", "d1 He","d2 He"),lty=c(1,1,2,2),col=c("black","red","black","red"),bty="n")
}

# plot sim results for many loci:
he.multi_plot<-function(d1sigma=NULL,d2sigma=NULL,Fphi=NULL,Mphi=NULL,he.data=NULL,N.data=NULL,nyears=tmax,...){

if (n>2) stop ("function not defined for > 2 demes")

# set data:
hedat<-he.data[,grep("rep",colnames(he.data))]
ndat<-N.data[,grep("rep",colnames(N.data))]

# make ylim offsets:
rangepopsize<-max(ndat)-min(ndat)
n.multiplier<-rangepopsize/max(ndat)
if(n.multiplier<0.2) noffset<-max(ndat)/6 else noffset<-max(ndat)/2

rangehe<-max(hedat)-min(hedat)
he.multiplier<-rangehe/max(hedat)
if(he.multiplier<0.2) heoffset<-max(hedat)/16 else heoffset<-max(hedat)/2

# set parameters and titles:
par(mar=c(4,4,3,4),mgp=c(2.5,1,0),new=F)
plot(1:nyears,1:nyears,type="n",col="gray80",xlab="Years",ylab="",bty="u",las=1,ylim=c(min(ndat)-noffset,max(ndat)),yaxt="n")
axis(side=2,at=axTicks(2),labels=axTicks(2)/1000,las=1)
title(main=bquote(sigma["d1"]~" = "~.(d1sigma)*"; "*sigma["d2"]~" = "~.(d2sigma)),line=2)
title(main=bquote(phi["F"]~" = "~.(Fphi)*"; "*phi["M"]~" = "~.(Mphi)),line=0)
title(ylab="Population size (,000 individuals)",mgp=c(2.5,1,0))

# plot population size:
for(i in 1:length(ndat)){

lines(1:nyears,ndat[1:nyears,i],col=rgb(0,0,0,0.15))
lines(1:nyears,ndat[(nyears+1):nrow(ndat),i],col=rgb(1,0,0,0.3))

} # close pop size

# reset y axis:
par(new=T)
plot(1:nyears,1:nyears,axes=F,type="n",xlab=NA,ylab="",ylim=c(min(hedat)-heoffset,max(hedat)+heoffset),las=1,lty=2)

# plot heterozygosity:
for(i in 1:length(hedat)){

lines(1:nyears,hedat[1:nyears,i],col=rgb(0,0,0,0.15),lty=2)
lines(1:nyears,hedat[(nyears+1):nrow(ndat),i],col=rgb(1,0,0,0.3),lty=2)

} # close he

# Add labels:
axis(side=4,las=1)
mtext("He",side=4,cex = par("cex.lab"),line=2.5,las=0)
legend("bottomleft",legend=c("d1 N","d2 N", "d1 He","d2 He"),lty=c(1,1,2,2),col=c("black","red","black","red"),bty="n")

} # close he.multi_plot

# plot sim results for many loci, putting multiple plots per page:
he.mplot_panel<-function(d1sigma=NULL,d2sigma=NULL,Fphi=NULL,Mphi=NULL,delta=NULL,he.data=NULL,N.data=NULL,nyears=tmax,...){

if (n>2) stop ("function not defined for > 2 demes")

# set data:
hedat<-he.data[,grep("rep",colnames(he.data))]
ndat<-N.data[,grep("rep",colnames(N.data))]

# make ylim offsets:
rangepopsize<-max(ndat)-min(ndat)
n.multiplier<-rangepopsize/max(ndat)
if(n.multiplier<0.2) noffset<-max(ndat)/6 else noffset<-max(ndat)/2

rangehe<-max(hedat)-min(hedat)
he.multiplier<-rangehe/max(hedat)
if(he.multiplier<0.2) heoffset<-max(hedat)/16 else heoffset<-max(hedat)/2

# set parameters and titles:
par(mar=c(4,4,3,4),mgp=c(2.5,1,0),new=F)
plot(1:nyears,1:nyears,type="n",col="gray80",xlab="Years",ylab="",bty="u",las=1,ylim=c(min(ndat)-noffset,max(ndat)),yaxt="n")
axis(side=2,at=axTicks(2),labels=axTicks(2)/1000,las=1)
title(main=bquote(sigma["d1"]~" = "~.(d1sigma)*"; "*sigma["d2"]~" = "~.(d2sigma)),line=1.8)
title(main=bquote(phi["F"]~" = "~.(Fphi)*"; "*phi["M"]~" = "~.(Mphi)*"; "*delta~" = "~.(delta)),line=0.5)
title(ylab="Population size (,000)",mgp=c(2.5,1,0))

# plot population size:
for(i in 1:length(ndat)){

lines(1:nyears,ndat[1:nyears,i],col=rgb(0,0,0,0.15))
lines(1:nyears,ndat[(nyears+1):nrow(ndat),i],col=rgb(1,0,0,0.3))

} # close pop size

# reset y axis:
par(new=T)
plot(1:nyears,1:nyears,axes=F,type="n",xlab=NA,ylab="",ylim=c(min(hedat)-heoffset,max(hedat)+heoffset),las=1,lty=2)

# plot heterozygosity:
for(i in 1:length(hedat)){

lines(1:nyears,hedat[1:nyears,i],col=rgb(0,0,0,0.15),lty=2)
lines(1:nyears,hedat[(nyears+1):nrow(ndat),i],col=rgb(1,0,0,0.3),lty=2)

} # close he

# Add labels:
axis(side=4,las=1,ylab="He")
mtext("He",cex=par("cex"),side=4,line=2.5,las=0)
# legend("bottomleft",legend=c("d1 N","d2 N", "d1 He","d2 He"),lty=c(1,1,2,2),col=c("black","red","black","red"),bty="n")

} # close he.multi_plot





