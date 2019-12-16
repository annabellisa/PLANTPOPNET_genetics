
# ------------------------------------ #
# ----------- SUPPLEMENT 06  --------- #
# ------------------------------------ #

### Summarise overlap among outlier detection methods:
### Author: Annabel Smith

# Load workspace:
load("../04_workspaces/Supp_venn_diagram")

# Load functions:
invisible(lapply(paste("../02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

library(VennDiagram)

# INSTRUCTIONS:

# run STEP_02_filter_and_format.R all the way down to the outlier loci step...

# remove cultivars and outgroups, remove sites with small sample sizes and run the script:

# Create data set with ONLY non-neutral loci
# Then save.image("../04_workspaces/Supp_venn_diagram")

# NEW version (after bioclim update): dim() should 454 x 3028 (3026 loci)

ghead(filtered_data); dim(filtered_data)

# Subset res for loci in filtered_data:
ladapt<-colnames(filtered_data)[3:length(filtered_data)]
head(ladapt); length(ladapt)

res_sub<-res[which(res$locus %in% ladapt),]
head(res_sub); dim(res_sub)

# Summarise overlap:
res2<-res_sub[,c("bs_outl","pca_outl","lfmm_outl")]
res2$sum1<-rowSums(cbind(res2[,1],res2[,2]))
res2$sum2<-rowSums(cbind(res2[,2],res2[,3]))
res2$sum3<-rowSums(cbind(res2[,1],res2[,3]))
res2$sum4<-rowSums(cbind(res2[,1],res2[,2],res2[,3]))
res2$sum1<-ifelse(res2$sum1==1,0,res2$sum1)
res2$sum2<-ifelse(res2$sum2==1,0,res2$sum2)
res2$sum3<-ifelse(res2$sum3==1,0,res2$sum3)
res2$sum4<-ifelse(res2$sum4==1,0,res2$sum4)
res2$sum4<-ifelse(res2$sum4==2,0,res2$sum4)
res2$sum1<-ifelse(res2$sum1==2,1,res2$sum1)
res2$sum2<-ifelse(res2$sum2==2,1,res2$sum2)
res2$sum3<-ifelse(res2$sum3==2,1,res2$sum3)
res2$sum4<-ifelse(res2$sum4==3,1,res2$sum4)
head(res2)

res3<-colSums(res2)
res3

venn.plot<-draw.triple.venn(res3[1],res3[2],res3[3],res3[4],res3[5],res3[6],res3[7], category=c("BayeScan","PCAdapt","LFMM"),fill=rgb(0,0,0,0.5),fontfamily="sans",cat.fontfamily="sans",cex=2)

quartz("xx",8,8,pointsize=18,dpi=80,file="xx.pdf",type="pdf")
grid.draw(venn.plot)
dev.off()







