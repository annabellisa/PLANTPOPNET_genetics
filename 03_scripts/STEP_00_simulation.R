
# ------------------------------------ #
# ------------- STEP 00  ------------- #
# ------------------------------------ #

### SIMULATION EXPERIMENTS

# Author: Annabel Smith

# Takes 52 mins on Big Mac with 300 reps

# Load functions (BIG MAC):
invisible(lapply(paste("/Volumes/LaCie/Annabel_GBS/sims_Feb_2019/sims_19Feb2019/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Load functions (laptop):
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Load workspace:
load("../04_workspaces/STEP00_sims.RData")

# Load packages:
library(MetaPopGen)
# install.packages("MetaPopGen.0.0.4.tar.gz",type="source",repos=NULL)

image.name<-paste("sims_",Sys.Date(),".RData", sep="")
save.image(image.name)

# ****
#   	PARAMETERISE	  # ----

# The model is not spatially explicit, but we modelled the number of individuals based on the density in an area of 1 ha (10,000 m2). 

# Jesus Villellas recommended using site means for density (mean = 103, range = 15:400 in PlantPopNet data) because the population doesn't operate at a scale as fine as the plot. 150,000 individuals total (25000 / age-class and genotype with age-structure data) is the lower end of density in the PPN data (15/m2):

nind <- 25000	# Initial number of individuals per age-class and genotype
n <- 2	# Number of demes
z <- 2	# Number of age classes
l <- 2	# Number of alleles
m <- 3	# Number of genotypes; this must be related to l as m<-l*(l+1)/2
tmax <- 100	# Number of time steps (generations) in the simulation

# Generation time in P. lanceolata is 2.78 (Steiner et al. 2018), or between 1 and 3 (van Groenendael & Slim, 1988). In the model, time steps are generations, so we adjust the number of time steps to achieve the relevant number of generations. We want to ensure converge of the estimates and also a timeframe relevant to the recent introduction history. We are basing the simulation on 300 years to account for the introduction history of non-native populations. We could expect P. lanceolata to have undergone 100 generations in that time. 

# Initialise variables, using the correct dimensions
sigma <- array(NA,dim=c(m,n,z,tmax))
phi_M <- array(0,dim=c(m,n,z,tmax))
phi_F <- array(0,dim=c(m,n,z,tmax))

# To get a reasonable value for survival, use a predicted value for longevity and the proportion of individuals estimated to be alive after that time period. E.g. 5 % of individuals will be alive after 5 years. van Groenendael & Slim (1998, J Ecol) report 5-10% alive after 5 years. 
longev<-5
prop_alive<-0.05
surv_rate<-round(exp(log(prop_alive)/longev),2)

# This also gives 5% alive after 5 years but assumes a seedling survival rate of 0.1 for the first year and an annual survival rate for adults over four years:

# d1.sig<-surv_rate
# d2.sig<-0.3
d1.sigj<-0.1 # i.e. seedlings
d1.siga<-(0.05/d1.sigj)^0.25 # i.e. adults
d2.sigj<-0.2 # i.e. seedlings
d2.siga<-(0.05/d2.sigj)^0.25 # i.e. adults

# check it:
0.1*0.84^4

# Values from PPN 2015 data. 
Fphi<-20
Mphi<-10000

# Assign values to variables
# sigma[,,,] <- d1.sig
sigma[,1,1,] <- d1.sigj
sigma[,1,2,] <- d1.siga
sigma[,2,1,] <- d2.sigj
sigma[,2,2,] <- d2.siga

# Give fecundities to the adult stage only, leaving the seedling stage at zero:
phi_F[,,2,] <- Fphi
phi_M[,,2,] <- Mphi

# Dispersal matrix
m.m <- 0.0 # no dispersal in Experiment 1
delta <- matrix(rep(m.m/(n-1),n^2),nrow=n,ncol=n)
diag(delta)<-rep((1-m.m),n)

# Mutation matrix
mu <- array(1e-6,dim=c(l,l))
mu[1,1] <- 1 - mu[2,1]
mu[2,2] <- 1 - mu[1,2]

# Maximum recruitment (nind*m will give the number of individuals in a single deme):
# kappa0 <- array(nind*m,c(n,tmax))

# Use Jesus's seedling survival rates (proportion of seeds becoming seedlings) to set the recruitment parameter. The total number of seedlings produced should be the female fecundity x population size, then multiply this by the germination rate:
g_rates<-c(0.0043,0.0049, 0.0206, 0.0082, 0.0229, 0.0341, 0.0851, 0.1292)
mean_grate<-mean(g_rates)
range(g_rates)
kappa0 <- array(round((Fphi*(nind*3))* mean_grate,0),c(n,tmax))

# close parameterise ---- 

start.time<-print(Sys.time())
save.image(image.name)

# ****
# EXPERIMENT 1: FEMALE FECUNDITY & SURVIVAL ----

# fem.fecund as below; no.reps = 300; initial pop size = 150,000 (15 / m2)

# Plot He over time as simulation runs:
dev.new(title="",width=11,height=8, dpi=80, pointsize=14, noRStudioGD = T,file=paste("exp1_time.pdf",sep=""), type="pdf")
par(mfrow=c(4,4))

fem.fecund<-c(1,5,10,15,20,25,30,35,40,45,50,70,90,100)

ff<-data.frame(fem.fecund=fem.fecund,deme=c(rep("d1",length(fem.fecund)),rep("d2",length(fem.fecund))),He=NA, He.se=NA, N=NA, N.se=NA)

# Effect of female fecundity on genetic diversity, no migration:

for (i in 1:length(fem.fecund)){
  
  # update extinct indicator:
  ext.ind<-data.frame(deme=c("d1","d2"),extinct=0)
  
  fecund.thisrun<-fem.fecund[i]
  
  # Update female fecundity:
  phi_F[,,2,] <- fecund.thisrun
  
  # Initial population structure:
  N1<-pop.wrapper1()
  
  # RUN SIMULATION:
  N<-sim.wrapper2(nrep=300)
  
  # Summarise heterozygosity and N:
  he.allruns<-he.multi(nrep=300,ndeme=n)
  N.allruns<-popsize.multi(nrep=300,ndeme=n)
  
  # If the population went extinct:
  if(length(which(is.na(he.allruns[he.allruns$deme=="d1",])))>0)  ext.ind$extinct[ext.ind$deme=="d1"]<-1
  if(length(which(is.na(he.allruns[he.allruns$deme=="d2",])))>0)  ext.ind$extinct[ext.ind$deme=="d2"]<-1
  
  # Summarise He and N at END of simulation:
  d1.he<-as.numeric(he.allruns[tail(which(he.allruns$deme=="d1"),1),-which(colnames(he.allruns)=="deme")])
  d2.he<-as.numeric(he.allruns[tail(which(he.allruns$deme=="d2"),1),-which(colnames(he.allruns)=="deme")])
  d1.N<-as.numeric(N.allruns[tail(which(N.allruns$deme=="d1"),1),-which(colnames(he.allruns)=="deme")])
  d2.N<-as.numeric(N.allruns[tail(which(N.allruns$deme=="d2"),1),-which(colnames(he.allruns)=="deme")])
  
  ff$He[ff$fem.fecund==fecund.thisrun & ff$deme=="d1"]<-mean(d1.he)
  ff$He[ff$fem.fecund==fecund.thisrun & ff$deme=="d2"]<-mean(d2.he)
  ff$He.se[ff$fem.fecund==fecund.thisrun & ff$deme=="d1"]<-std.err(d1.he)
  ff$He.se[ff$fem.fecund==fecund.thisrun & ff$deme=="d2"]<-std.err(d2.he)
  ff$N[ff$fem.fecund==fecund.thisrun & ff$deme=="d1"]<-mean(d1.N)
  ff$N[ff$fem.fecund==fecund.thisrun & ff$deme=="d2"]<-mean(d2.N)
  ff$N.se[ff$fem.fecund==fecund.thisrun & ff$deme=="d1"]<-std.err(d1.N)
  ff$N.se[ff$fem.fecund==fecund.thisrun & ff$deme=="d2"]<-std.err(d2.N)
  
  # Plot results over time:
  
  # If either population went extinct, make all NAs zero for plotting:
  if(length(which(ext.ind$extinct==1))>0) he.allruns[is.na(he.allruns)]<-0
  
  he.mplot_panel(d1sigma=paste(sigma[1,1,1,1],round(sigma[1,1,2,1],2),sep="/"),d2sigma=paste(sigma[1,2,1,1],round(sigma[1,2,2,1],2),sep="/"),Fphi=phi_F[1,1,2,1],Mphi=phi_M[1,1,2,1],delta=0,he.data=he.allruns,N.data=N.allruns,nyears=tmax)
  
  if(ext.ind$extinct[1]==1) text(50,0.4,"d1 extinct",adj=0,cex=1.5)
  if(ext.ind$extinct[2]==1) text(50,0.2,"d2 extinct",adj=0,cex=1.5)

} # close for i

dev.off()

# Add CIs to data:

ff$He.lci<-ff$He-(1.96*ff$He.se)
ff$He.uci<-ff$He+(1.96*ff$He.se)
ff$N.lci<-ff$N-(1.96*ff$N.se)
ff$N.uci<-ff$N+(1.96*ff$N.se)
# write.table(ff,"ff.txt",sep="\t",row.names=F,quote=F)

# PLOT RESULTS:

dev.new(title="",width=13,height=6,pointsize=22,dpi=80,noRStudioGD = T,file=paste("exp1_result.pdf",sep=""), type="pdf")
# dev.new(title="",width=13,height=6,pointsize=22,dpi=80,noRStudioGD = T)
par(mfrow=c(1,2),mar=c(4,4.5,0.5,2),mgp=c(3.5,1,0))

plot(ff$fem.fecund[ff$deme=="d1"],ff$He[ff$deme=="d1"],type="l",ylim=c(min(ff$He.lci,na.rm=T),max(ff$He.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="Heterozygosity",lwd=2)
title(xlab="Female fecundity (seeds / plant)",mgp=c(2.5,1,0))
xx<-ff$fem.fecund[ff$deme=="d1"]
lci<-ff$He.lci[ff$deme=="d1"]
uci<-ff$He.uci[ff$deme=="d1"]
xx<-xx[-which(is.na(lci))]
lci<-lci[-which(is.na(lci))]
uci<-uci[-which(is.na(uci))]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,0,0.15), border=NA)

lines(ff$fem.fecund[ff$deme=="d2"],ff$He[ff$deme=="d2"],col="red",lwd=2)
xx<-ff$fem.fecund[ff$deme=="d2"]
lci<-ff$He.lci[ff$deme=="d2"]
uci<-ff$He.uci[ff$deme=="d2"]
xx<-xx[-which(is.na(lci))]
lci<-lci[-which(is.na(lci))]
uci<-uci[-which(is.na(uci))]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(1,0,0,0.15), border=NA)

text(ff$fem.fecund[which(is.na(ff$He[ff$deme=="d1"]))],0.465,"†",cex=1.1)
text(ff$fem.fecund[which(is.na(ff$He[ff$deme=="d2"]))],0.49,"†",cex=1.1,col="red")

legend("bottomright",legend=c(as.expression(bquote(sigma["d1"]~" = "~.(paste(sigma[1,1,1,1],round(sigma[1,1,2,1],2),sep="/")))),bquote(sigma["d2"]~" = "~.(paste(sigma[1,2,1,1],round(sigma[1,2,2,1],2),sep="/"))),"  †   =  extinct"),col=c("black","red","white"),lty=1,bty="n",lwd=2)

plot(ff$fem.fecund[ff$deme=="d1"],ff$N[ff$deme=="d1"],type="l",ylim=c(min(ff$N.lci,na.rm=T),max(ff$N.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="",yaxt="n",lwd=2)
axis(side=2,at=axTicks(2),labels=axTicks(2)/1000,las=1)
title(ylab="Population size (,000)",mgp=c(2.5,1,0))
title(xlab="Female fecundity (seeds / plant)",mgp=c(2.5,1,0))
xx<-ff$fem.fecund[ff$deme=="d1"]
lci<-ff$N.lci[ff$deme=="d1"]
uci<-ff$N.uci[ff$deme=="d1"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,0,0.15), border=NA)

lines(ff$fem.fecund[ff$deme=="d2"],ff$N[ff$deme=="d2"],col="red",lwd=2)
xx<-ff$fem.fecund[ff$deme=="d2"]
lci<-ff$N.lci[ff$deme=="d2"]
uci<-ff$N.uci[ff$deme=="d2"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(1,0,0,0.15), border=NA)

legend("bottomright",legend=c("Deme 1", "Deme 2"),col=c("black","red"),lty=1,bty="n",lwd=2)

dev.off()

# mtext("Density = 15 / m2, no migration, 300 year simulation, 30 loci",side=3,adj=0,font=1,at=-100,line=1.5)

save.image(image.name)

# close EXP1 ----

# ****
# EXPERIMENT 2: MIGRATION ----

# from the Andrello paper, the migration parameter is: the probabilities that a newborn in deme r will successfully disperse to deme i at time t

# Dispersal happens before recruitment, so this is the total number of new borns. In four of five age classes, each of 25000 plants produce 20 seeds. This gives 500,000 newborns per genotype and 1,500,000 per deme. The kappa parameter specifies that only 57,000 can recruit but this happens after dispersal. 

# So if we have 1,500,000 newborns and the migration rate is 1%, we still have 15,000 individuals migrating, which is a huge number.

newborns<-1500000
migration_rate<-0.01
mig_percent<-migration_rate*100
mig_percent

newborns*migration_rate

# 0.0001 % migration rate gives 1.5 per generation
# 0.00001 % migration rate gives 0.15 per generation

# Initially ran the range from 0.0000001 (0.15/generation) to 0.00001 (15/generation), 0.0001 (150/generation) and 0.001 (1500/generation) and it DID NOT reach the mixing threshold

# Then I ran 0.0000001 (0.15/generation) to 0.01 (15000/generation)...

# fem.fecund = 20; no.reps = 30; initial pop size = 150,000 (15 / m2); mig_rates<-c(0,seq(0.0000001,0.001,length.out=10))

mig_rates*newborns

# Re-set parameters that were changed in the previous experiment:
phi_F[,,2,] <- Fphi

dev.new(title="",width=11,height=8, dpi=80, pointsize=14, noRStudioGD = T,file=paste("exp2_time.pdf",sep=""), type="pdf")
par(mfrow=c(4,4))

# setting it by the number of migrants per generation and then dividing by newborns to get the migration rate:
mig_rates<-c(0,1,exp(1:11))/newborns

mr<-data.frame(mig.rate=mig_rates,deme=c(rep("d1",length(mig_rates)),rep("d2",length(mig_rates))),He=NA, He.se=NA, N=NA, N.se=NA)

mig.diff<-data.frame(mig.rate=mig_rates,He_diff=NA,diff_se=NA)

# Effect of migration on difference in genetic diversity:

for (i in 1:length(mig_rates)){
  
  # update extinct indicator:
  ext.ind<-data.frame(deme=c("d1","d2"),extinct=0)
  
  # update migration rate:
  mig.thisrun<-mig_rates[i]
  
  m.m <- mig.thisrun
  delta <- matrix(rep(m.m/(n-1),n^2),nrow=n,ncol=n)
  diag(delta)<-rep((1-m.m),n)
  
  # Initial population structure:
  N1<-pop.wrapper1()
  
  # RUN SIMULATION:
  N<-sim.wrapper2(nrep=300)
  
  # Summarise heterozygosity and N:
  he.allruns<-he.multi(nrep=300,ndeme=n)
  N.allruns<-popsize.multi(nrep=300,ndeme=n)
  
  # If the population went extinct:
  if(length(which(is.na(he.allruns[he.allruns$deme=="d1",])))>0)  ext.ind$extinct[ext.ind$deme=="d1"]<-1
  if(length(which(is.na(he.allruns[he.allruns$deme=="d2",])))>0)  ext.ind$extinct[ext.ind$deme=="d2"]<-1
  
  # Summarise He and N at END of simulation:
  d1.he <-as.numeric(he.allruns[tail(which(he.allruns$deme=="d1"),1),-which(colnames(he.allruns)=="deme")])
  d2.he<-as.numeric(he.allruns[tail(which(he.allruns$deme=="d2"),1),-which(colnames(he.allruns)=="deme")])
  d1.N<-as.numeric(N.allruns[tail(which(N.allruns$deme=="d1"),1),-which(colnames(he.allruns)=="deme")])
  d2.N<-as.numeric(N.allruns[tail(which(N.allruns$deme=="d2"),1),-which(colnames(he.allruns)=="deme")])
  
  mig.diff$He_diff[i]<-mean(d1.he-d2.he)
  mig.diff$diff_se[i]<-std.err(d1.he-d2.he)
  
  mr$He[mr$mig.rate==mig.thisrun & mr$deme=="d1"]<-mean(d1.he)
  mr$He[mr$mig.rate==mig.thisrun & mr$deme=="d2"]<-mean(d2.he)
  mr$He.se[mr$mig.rate==mig.thisrun & mr$deme=="d1"]<-std.err(d1.he)
  mr$He.se[mr$mig.rate==mig.thisrun & mr$deme=="d2"]<-std.err(d2.he)
  mr$N[mr$mig.rate==mig.thisrun & mr$deme=="d1"]<-mean(d1.N)
  mr$N[mr$mig.rate==mig.thisrun & mr$deme=="d2"]<-mean(d2.N)
  mr$N.se[mr$mig.rate==mig.thisrun & mr$deme=="d1"]<-std.err(d1.N)
  mr$N.se[mr$mig.rate==mig.thisrun & mr$deme=="d2"]<-std.err(d2.N)
  
  # Plot results over time:
  
  # If either population went extinct, make all NAs zero for plotting:
  if(length(which(ext.ind$extinct==1))>0) he.allruns[is.na(he.allruns)]<-0
  
  he.mplot_panel(d1sigma=paste(sigma[1,1,1,1],round(sigma[1,1,2,1],2),sep="/"),d2sigma=paste(sigma[1,2,1,1],round(sigma[1,2,2,1],2),sep="/"),Fphi=phi_F[1,1,2,1],Mphi=phi_M[1,1,2,1],delta=round(mig.thisrun,6),he.data=he.allruns,N.data=N.allruns,nyears=tmax)
  
  if(ext.ind$extinct[1]==1) text(50,0.4,"d1 extinct",adj=0,cex=1.5)
  if(ext.ind$extinct[2]==1) text(50,0.2,"d2 extinct",adj=0,cex=1.5)

} # close i for

dev.off()

# Add CIs to data:

mr$He.lci<-mr$He-(1.96*mr$He.se)
mr$He.uci<-mr$He+(1.96*mr$He.se)
mr$N.lci<-mr$N-(1.96*mr$N.se)
mr$N.uci<-mr$N+(1.96*mr$N.se)

mig.diff$diff.lci<-mig.diff$He_diff-(1.96* mig.diff$diff_se)
mig.diff$diff.uci<-mig.diff$He_diff+(1.96* mig.diff$diff_se)

# Add number of migrants to mr and mig.diff:
mr$no_mig<-mr$mig.rate*newborns
mig.diff$no_mig<-mig.diff$mig.rate*newborns

# write.table(mr,"mr.txt",sep="\t",row.names=F,quote=F)
# write.table(mig.diff,"mig.txt",sep="\t",row.names=F,quote=F)

# PLOT RESULTS:

dev.new(title="",width=13,height=6,pointsize=22,dpi=80,noRStudioGD = T,file=paste("exp2_result.pdf",sep=""), type="pdf")
# dev.new(title="",width=13,height=6,pointsize=22,dpi=80,noRStudioGD = T)
par(mfrow=c(1,2),mar=c(4,4.5,0.5,2),mgp=c(3.5,1,0))

xl<-round(mr$no_mig[mr$deme=="d1"],0)
xl<-c(xl[1:5],round(xl[6:8],-1),round(xl[9:13],-3))

plot(log(mr$no_mig[mr$deme=="d1"]),mr$He[mr$deme=="d1"],type="l",ylim=c(min(mr$He.lci,na.rm=T),max(mr$He.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="Heterozygosity",xaxt="n",lwd=2)
axis(side=1,at=log(mr$no_mig[mr$deme=="d1"]),labels=xl)
title(xlab="Migrants / generation (log scale)",mgp=c(2.5,1,0))
xx<-log(mr$no_mig[mr$deme=="d1"])
lci<-mr$He.lci[mr$deme=="d1"]
uci<-mr$He.uci[mr$deme=="d1"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,0,0.15), border=NA)

lines(log(mr$no_mig[mr$deme=="d2"]),mr$He[mr$deme=="d2"],col="red",lwd=2)
xx<-log(mr$no_mig[mr$deme=="d2"])
lci<-mr$He.lci[mr$deme=="d2"]
uci<-mr$He.uci[mr$deme=="d2"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(1,0,0,0.15), border=NA)

legend("bottomright",legend=c(as.expression(bquote(sigma["d1"]~" = "~.(paste(sigma[1,1,1,1],round(sigma[1,1,2,1],2),sep="/")))),bquote(sigma["d2"]~" = "~.(paste(sigma[1,2,1,1],round(sigma[1,2,2,1],2),sep="/")))),col=c("black","red"),lty=1,bty="n",lwd=2)

plot(log(mig.diff$no_mig),mig.diff$He_diff,pch=20,ylim=c(min(mig.diff$diff.lci,na.rm=T),max(mig.diff$diff.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="",xaxt="n")
axis(side=1,at=log(mig.diff$no_mig),labels=xl)
title(ylab="Difference in heterozygosity",mgp=c(3.5,1,0))
title(xlab="Migrants / generation (log scale)",mgp=c(2.5,1,0))
arrows(log(mig.diff$no_mig),mig.diff$diff.lci,log(mig.diff$no_mig),mig.diff$diff.uci,length=0.02,code=3,angle=90)
arrows(-15,0,log(max(mig.diff$no_mig)*2),0,length=0.02,code=3,angle=90)

dev.off()

save.image(image.name)

end.time<-print(Sys.time())
start.time
end.time
 
# close EXP2 ----

# ****
# PLOTS FOR PAPER 

dev.new(title="",width=10,height=10,pointsize=22,dpi=70,noRStudioGD = T)
par(mfrow=c(2,2),mar=c(4,5,2,1),mgp=c(3.5,1,0))

# EXPERIMENT 1: PLOT ----

plot(ff$fem.fecund[ff$deme=="d1"],ff$He[ff$deme=="d1"],type="l",ylim=c(min(ff$He.lci,na.rm=T),max(ff$He.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="Heterozygosity",lwd=2)
title(xlab="Female fecundity (seeds / plant)",mgp=c(2.2,1,0))
xx<-ff$fem.fecund[ff$deme=="d1"]
lci<-ff$He.lci[ff$deme=="d1"]
uci<-ff$He.uci[ff$deme=="d1"]
xx<-xx[-which(is.na(lci))]
lci<-lci[-which(is.na(lci))]
uci<-uci[-which(is.na(uci))]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,0,0.15), border=NA)

lines(ff$fem.fecund[ff$deme=="d2"],ff$He[ff$deme=="d2"],col="blue",lwd=2,lty=2)
xx<-ff$fem.fecund[ff$deme=="d2"]
lci<-ff$He.lci[ff$deme=="d2"]
uci<-ff$He.uci[ff$deme=="d2"]
xx<-xx[-which(is.na(lci))]
lci<-lci[-which(is.na(lci))]
uci<-uci[-which(is.na(uci))]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,1,0.15), border=NA)

mtext("(a)",side=3,at=-10,line=0.5)

text(ff$fem.fecund[which(is.na(ff$He[ff$deme=="d1"]))],0.487,"†",cex=1.1)
text(ff$fem.fecund[which(is.na(ff$He[ff$deme=="d2"]))],0.49,"†",cex=1.1,col="blue")

legend(20,0.4914,legend=c("juvenile survival:"),lty=c(1),col=rgb(0,0,0,0),bty="n",lwd=2, cex=0.9)
legend(41,0.49,legend=c("low","high"),col=c("black","blue"),lty=c(1,2),bty="n",lwd=2, cex=0.9)
legend(54.2,0.4874,legend=c("extinct"),col=c("black"),pch="†",bty="n", cex=0.9)

plot(ff$fem.fecund[ff$deme=="d1"],ff$N[ff$deme=="d1"],type="l",ylim=c(min(ff$N.lci,na.rm=T),max(ff$N.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="",yaxt="n",lwd=2)
axis(side=2,at=axTicks(2),labels=axTicks(2)/1000,las=1)
title(ylab="Population size (,000)",mgp=c(2.5,1,0))
title(xlab="Female fecundity (seeds / plant)",mgp=c(2.2,1,0))
xx<-ff$fem.fecund[ff$deme=="d1"]
lci<-ff$N.lci[ff$deme=="d1"]
uci<-ff$N.uci[ff$deme=="d1"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,0,0.15), border=NA)

lines(ff$fem.fecund[ff$deme=="d2"],ff$N[ff$deme=="d2"],col="blue",lwd=2, lty=2)
xx<-ff$fem.fecund[ff$deme=="d2"]
lci<-ff$N.lci[ff$deme=="d2"]
uci<-ff$N.uci[ff$deme=="d2"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(1,0,0,0.15), border=NA)

mtext("(b)",side=3,at=-10,line=0.5)

legend(20,24000,legend=c("juvenile survival:"),lty=c(1),col=rgb(0,0,0,0),bty="n",lwd=2, cex=0.9)
legend(41,18500,legend=c("low","high"),col=c("black","blue"),lty=c(1,2),bty="n",lwd=2, cex=0.9)

# close exp1 plot ----

# EXPERIMENT 2: PLOT ----

xl<-round(mr$no_mig[mr$deme=="d1"],0)
xl<-c(xl[1:5],round(xl[6:8],-1),round(xl[9:13],-3))

plot(log(mr$no_mig[mr$deme=="d1"]),mr$He[mr$deme=="d1"],type="l",ylim=c(0.496,max(mr$He.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="Heterozygosity",yaxt="n",xaxt="n",lwd=2)
axis(side=1,at=log(mr$no_mig[mr$deme=="d1"]),labels=xl)
axis(side=2, at=c(0.496,0.497,0.498,0.499),las=1)
title(xlab="Migrants / generation (log scale)",mgp=c(2.2,1,0))
xx<-log(mr$no_mig[mr$deme=="d1"])
lci<-mr$He.lci[mr$deme=="d1"]
uci<-mr$He.uci[mr$deme=="d1"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,0,0.15), border=NA)

lines(log(mr$no_mig[mr$deme=="d2"]),mr$He[mr$deme=="d2"],col="blue",lwd=2,lty=2)
xx<-log(mr$no_mig[mr$deme=="d2"])
lci<-mr$He.lci[mr$deme=="d2"]
uci<-mr$He.uci[mr$deme=="d2"]
xvec <- c(xx, tail(xx, 1), rev(xx), xx[1])
yvec <- c(lci, tail(uci, 1), rev(uci), lci[1])
polygon(xvec, yvec, col=rgb(0,0,1,0.15), border=NA)

mtext("(c)",side=3,at=-1.15,line=0.5)

par(xpd=NA)
legend(2.7,0.49709,legend=c("juvenile survival:"),lty=c(1),col=rgb(0,0,0,0),bty="n",lwd=2, cex=0.9)
legend(5,0.4968,legend=c("low","high"),col=c("black","blue"),lty=c(1,2),bty="n",lwd=2, cex=0.9)
par(xpd=F)

plot(log(mig.diff$no_mig),mig.diff$He_diff,pch=20,ylim=c(min(mig.diff$diff.lci,na.rm=T),max(mig.diff$diff.uci,na.rm=T)),bty="l",las=1,xlab="",ylab="",xaxt="n")
axis(side=1,at=log(mig.diff$no_mig),labels=xl)
title(ylab="Difference in heterozygosity",mgp=c(4,1,0))
title(xlab="Migrants / generation (log scale)",mgp=c(2.2,1,0))
arrows(log(mig.diff$no_mig),mig.diff$diff.lci,log(mig.diff$no_mig),mig.diff$diff.uci,length=0.02,code=3,angle=90)
arrows(-15,0,log(max(mig.diff$no_mig)*2),0,length=0.02,code=3,angle=90)

mtext("(d)",side=3,at=-1.15,line=0.5)


#  close exp1 plot ----






















