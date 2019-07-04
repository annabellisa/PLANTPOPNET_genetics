# Analyse FST for outgroups, then REMOVE:
### Author: Annabel Smith

gp_dir<-"/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSIS_RESULTS/Genepop_DATA_FILES"
dir(gp_dir)

ds_filt2<-read.table(paste(gp_dir,"dartseq_filt2.txt",sep="/"),header=T)
ghead(ds_filt2); dim(ds_filt2)

snp<-ds_filt2[,3:length(ds_filt2)]
snp<-data.frame(t(snp))
colnames(snp)<-ds_filt2$ind
ghead(snp)

typed_at<-apply(snp,2,function(x) 1-length(which(is.na(x)))/nrow(snp))
typed_at<-data.frame(site=ds_filt2$site,ind=names(typed_at),typed_at=typed_at)
typed_at<-typed_at[order(typed_at$ind),]
typed_at<-tidy.df(typed_at)
head(typed_at)

site_typed<-data.frame(site=names(tapply(typed_at$typed_at,typed_at$site,mean)),mean_typed=tapply(typed_at$typed_at,typed_at$site,mean))
site_typed<-tidy.df(site_typed)
# plot(site_typed,cex.axis=0.4,las=2)
head(site_typed)
snp<-NULL

mean_OG<-mean(site_typed$mean_typed[grep("OG",site_typed$site)])
mean_lanceolata<-mean(site_typed$mean_typed[-grep("OG",site_typed$site)])
save.image("../04_workspaces/STEP05_env_wksp")

head(pw,2)

og_lines<-c(grep("OG",pw$pop1),grep("OG",pw$pop2))
ogc_lines<-c(grep("OGC",pw$pop1),grep("OGC",pw$pop2))
ogm_lines<-c(grep("OGM",pw$pop1),grep("OGM",pw$pop2))

# quartz("",4,6)
# par(mfrow=c(3,2),mgp=c(2.5,1,0),mar=c(4,4,2,1))
# plot(1:nrow(site_typed),site_typed$mean_typed,cex.axis=1,las=2,xlab="",xaxt="n",ylab="Propn. genotyped",pch=20,col=ifelse(rownames(site_typed) %in% grep("OG",site_typed$site),"red","black"),main="Call rate",font.main=1,cex.main=1)
# legend(-1,0.45,legend=c("lanceolata","outgroup"),pch=20,col=c("black","red"),bty="n")
# hist(pw$fst[-og_lines],xlim=c(0,1),main=expression(italic("Plantago lanceolata")),xlab="Pairwise FST",ylab="",cex.main=1,font.main=1)
# hist(pw$fst[og_lines],xlim=c(0,1),main="All outgroups",xlab="Pairwise FST",ylab="",cex.main=1,font.main=1)
# hist(pw$fst[ogm_lines],xlim=c(0,1),main=expression(italic("Plantago major")),xlab="Pairwise FST",ylab="",cex.main=1,font.main=1)
# text(0.85,15,"vs. P. lanceolata",srt=90,adj=0)
# text(0.05,15,"vs. other P. major",srt=0,adj=0)

# hist(pw$fst[ogc_lines],xlim=c(0,1),main=expression(italic("Plantago coronopus")),xlab="Pairwise FST",ylab="",cex.main=1,font.main=1)
# text(0.95,10,"vs. P. major",srt=90,adj=0)
# text(0.55,10,"vs. P. lanceolata",srt=90,adj=0)
# boxplot(pw$fst[-og_lines],pw$fst[ogm_lines],pw$fst[ogc_lines],ylim=c(0,1),las=1,ylab="Pairwise FST")
# axis(side=1,at=1:3,labels=c("lan","maj","cor"))



