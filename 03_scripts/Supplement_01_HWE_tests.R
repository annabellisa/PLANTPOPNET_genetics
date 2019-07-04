
# ------------------------------------ #
# ----------- SUPPLEMENT 02  --------- #
# ------------------------------------ #

### HWE TESTS
### Author: Annabel Smith

# Load formatted data including including functions and c("snp_onerow","linf","sdat","ukril_snp")
load("../04_workspaces/STEP02_filter_wksp")

# load functions:
invisible(lapply(paste("/Users/annabelsmith/Documents/01_Current/PROJECTS/01_PLANTPOPNET/DATA_and_ANALYSIS/SNP_analysis/GENOTYPE_processing/ANALYSE_GENOME/02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Subset UK/IRL data, which form a single genetic cluster, for LD testing
ukirl_sites<-as.character(sdat$site_code[sdat$country %in% c("UK","Ireland")])
ukirl_sites<-ukirl_sites[-which(ukirl_sites %in% c("OGM5","CCT"))]
ukirl_snp<-snp_onerow[snp_onerow$site %in% ukirl_sites,]
ukirl_snp<-tidy.df(ukirl_snp)
ghead(ukirl_snp)
levels(ukirl_snp$site)

## --- *** Run HWE tests *** --- ##

# First, run the HWE test on the UK/IRL sample:
uki_filt1<-mono_loci(ukirl_snp,3)
uki_hwe<-hwe_exact(uki_filt1)
# write.table(uki_hwe,"uki_hwe.txt",quote=F,row.names=F,sep="\t")

# 29000 loci do not exist in this data, set, so we need some way to look for deviations in the overall sample. 

# Remove cultivars, outgroups and mono:
filt_ogcult<-snp_onerow[-c(grep("OG",snp_onerow$site),grep("CAT",snp_onerow$site),grep("CCT",snp_onerow$site),grep("CTP",snp_onerow$site)),]
filt_ogcult<-tidy.df(filt_ogcult)

# Remove any site with 7 or fewer individual
filt_sml<-filt_ogcult[-which(filt_ogcult$site %in% names(table(all_nomono$site)[table(all_nomono$site)<8])),]
filt_sml<-tidy.df(filt_sml)
ghead(filt_sml)
all_nomono<-mono_loci(filt_sml,3)
ghead(all_nomono)

# Then run the test for each population separately:

sites_to_test<-levels(all_nomono$site)
out.list<-list()

for (i in 1:length(sites_to_test)){

site.thisrun<-sites_to_test[i]
data.thisrun<-all_nomono[all_nomono$site==site.thisrun,]
data.thisrun<-tidy.df(data.thisrun)
data.thisrun<-mono_loci(data.thisrun,3)
ghead(data.thisrun)
hwe.thisrun<-hwe_exact(data.thisrun)

out.list[[i]]<-data.frame(site=site.thisrun,loci_out=hwe.thisrun$locus[hwe.thisrun$p<0.05])

} # close for

hwe_res<-do.call(rbind, out.list)
head(hwe_res)







