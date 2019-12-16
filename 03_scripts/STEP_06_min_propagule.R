
# ------------------------------------ #
# ------------- STEP 06  ------------- #
# ------------------------------------ #

### Minimum number of propagules required to produce all genetic diversity in non-native regions

### Author: Annabel Smith

load("../04_workspaces/STEP06_propmin_wksp")

# load functions:
invisible(lapply(paste("../02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

# Data:
head(site); dim(site)
head(assig53); dim(assig53)
head(all53); dim(all53)

#########################################
####  	  RESULT SUMMARY:			 ####
#########################################

# 1000 calculations, WITH replacement (i.e. individuals could have been sampled twice):
region_private
head(propagule_result); dim(propagule_result)

# 10000 calculations, WITH replacement (i.e. individuals could have been sampled twice):
region_private2
head(propagule_result2); dim(propagule_result2)

# 1000 calculations, WITHOUT replacemtent (i.e. individuals can NOT have been sampled twice):
region_private_norepl
head(propagule_result3); dim(propagule_result3)

# Does replacement approach affect the result?
# Only a slight difference for the NZ data set, which might be expected for such a small sample size (9)
t.test(propagule_result[,2],propagule_result3[,2])
t.test(propagule_result[,3],propagule_result3[,3])
t.test(propagule_result[,4],propagule_result3[,4])
t.test(propagule_result[,5],propagule_result3[,5])
t.test(propagule_result[,6],propagule_result3[,6])

# Does number of reps affect the result?
# NO, so reporting 1000 is fine
t.test(propagule_result[,2],propagule_result2[,2])
t.test(propagule_result[,3],propagule_result2[,3])
t.test(propagule_result[,4],propagule_result2[,4])
t.test(propagule_result[,5],propagule_result2[,5])
t.test(propagule_result[,6],propagule_result2[,6])


#########################################
####  		  DATA SET-UP:			 ####
#########################################

# Make separate allele data frames:

# The genetic data are the same as used in all other neutral analyses:
filt_dat<-filtered_data[,3:length(filtered_data)]
ghead(filt_dat)

# Single row data:
# 0 = Reference allele homozygote (0101)
# 1 = SNP allele homozygote (0202)
# 2 = heterozygote (0102)

# Does it have the allele (1) or not (0)
ref_allele<-convert_geno(filt_dat, 2, 1, 0, 1, 0, 1, NA)
snp_allele<-convert_geno(filt_dat, 2, 1, 0, 1, 1, 0, NA)
colnames(ref_allele)<-paste(colnames(ref_allele),"_a", sep="")
colnames(snp_allele)<-paste(colnames(snp_allele),"_b", sep="")
ghead(ref_allele); dim(ref_allele)
ghead(snp_allele); dim(snp_allele)

both_alleles<-cbind(filtered_data[,1:2],ref_allele, snp_allele)
ghead(both_alleles); dim(both_alleles)

# Continental analysis for non-native range:

head(site); dim(site)
levels(site$region)

# sites for each region:
ntham_sites<-as.character(site$site_code[site$region=="Nth_America"])
afr_sites<-as.character(site$site_code[site$region=="Africa"])
nz_sites<-as.character(site$site_code[site$country=="NewZealand"])
jp_sites<-as.character(site$site_code[site$country=="Japan"])
au_sites<-as.character(site$site_code[site$country=="Australia"])

# snp data for each region:
ntham<-both_alleles[which(both_alleles$site %in% ntham_sites),]
ntham<-tidy.df(ntham)
afr<-both_alleles[which(both_alleles$site %in% afr_sites),]
afr<-tidy.df(afr)
nz<-both_alleles[which(both_alleles$site %in% nz_sites),]
nz<-tidy.df(nz)
jp<-both_alleles[which(both_alleles$site %in% jp_sites),]
jp<-tidy.df(jp)
au<-both_alleles[which(both_alleles$site %in% au_sites),]
au<-tidy.df(au)

# Assignment probabilities for the native range:
nat_dat<-all53[all53$native=="native",]
nat_dat<-tidy.df(nat_dat)
head(nat_dat); dim(nat_dat)

# SNP data for all of EUR:
eur_indivs<-as.character(nat_dat$indiv)
eur_snp<-both_alleles[which(both_alleles$ind %in% eur_indivs),]
eur_snp<-tidy.df(eur_snp)
ghead(eur_snp); dim(eur_snp)

# Alleles for each non-native region:

# This function returns a character vector of all alleles present in a specified snp dataframe:

get_alleles<-function(snp_dat){

alleles.thisrun<-data.frame(allele=names(colSums(snp_dat[,3:length(snp_dat)], na.rm=T)), allele_present=ifelse(as.numeric(colSums(snp_dat[,3:length(snp_dat)], na.rm=T))==0,0,1))
al<-as.character(alleles.thisrun$allele[alleles.thisrun$allele_present==1])
return(al)

} # close allele function

ghead(ntham); dim(ntham)
ghead(afr); dim(afr)
ghead(nz); dim(nz)
ghead(jp); dim(jp)
ghead(au); dim(au)

ntham_al<-get_alleles(ntham)
afr_al<-get_alleles(afr)
nz_al<-get_alleles(nz)
jp_al<-get_alleles(jp)
au_al<-get_alleles(au)
head(ntham_al); length(ntham_al)
head(afr_al); length(afr_al)
head(nz_al); length(nz_al)
head(jp_al); length(jp_al)
head(au_al); length(au_al)

# And for EUR:
eur_al<-get_alleles(eur_snp)
head(eur_al); length(eur_al)

#########################################
####  	  PRIVATE ALLELES:			 ####
#########################################

# Summarise unique alleles in each non-native region, relative to the whole EUR data set:
ntham_unique_EUR<-ntham_al[which(!ntham_al %in% eur_al)]
afr_unique_EUR<-afr_al[which(!afr_al %in% eur_al)]
nz_unique_EUR<-nz_al[which(!nz_al %in% eur_al)]
jp_unique_EUR<-jp_al[which(!jp_al %in% eur_al)]
au_unique_EUR<-au_al[which(!au_al %in% eur_al)]

region_private<-data.frame(region=c("NorthAmerica","Africa","NewZealand", "Japan", "Australia"), alleles_not_in_EUR=c(length(ntham_unique), length(afr_unique), length(nz_unique), length(jp_unique), length(au_unique)))

# Get private alleles for all populations:

ghead(both_alleles); dim(both_alleles)
head(site); dim(site)

# All alleles == n.loci*2 == 17160*2
all_al<-get_alleles(both_alleles)
head(all_al); length(all_al)

priv_al<-data.frame(site_code=site$site_code, private_alleles=NA)

for (i in 1:nrow(priv_al)){

site.thisrun<-as.character(priv_al$site_code[i])
snp.thisrun<-both_alleles[both_alleles$site==site.thisrun,]
snp.thisrun<-tidy.df(snp.thisrun)
ghead(snp.thisrun); dim(snp.thisrun)

al.thisrun<-get_alleles(snp.thisrun)
head(al.thisrun); length(al.thisrun)

snp.othersites<-both_alleles[-which(both_alleles$site==site.thisrun),]
al.othersites<-get_alleles(snp.othersites)
head(al.othersites); length(al.othersites)

unique.thisrun<-al.thisrun[which(!al.thisrun %in% al.othersites)]

priv_al$private_alleles[i]<-length(unique.thisrun)

} # close i

save.image("../04_workspaces/STEP06_propmin_wksp")

# write.table(priv_al, file="priv.txt", sep="\t", row.names=F, quote=F)

#########################################
####  	 MIN. NO. PROPAGULES:		 ####
#########################################

# Non-native alleles that are represented in the native range (use these for finding min. no. propagules, there's no point in including the unique alleles because it will never reach coverage):

ntham_rep<-ntham_al[which(ntham_al %in% eur_al)]
afr_rep<-afr_al[which(afr_al %in% eur_al)]
nz_rep<-nz_al[which(nz_al %in% eur_al)]
jp_rep<-jp_al[which(jp_al %in% eur_al)]
au_rep<-au_al[which(au_al %in% eur_al)]

# Total number of alleles:
region_private$total_no_alleles<-c(length(ntham_al),length(afr_al),length(nz_al),length(jp_al),length(au_al))
region_private$rep_alleles<-c(length(ntham_rep),length(afr_rep),length(nz_rep),length(jp_rep),length(au_rep))

# check: this should == the alleles not in EUR:
region_private$total_no_alleles-region_private$rep_alleles

# How many individuals from all of Europe would it take to represent all alleles in each of the non-native regions, disregarding the alleles that don't exist in Europe?

# These are the non-native alleles that can be represented by Europe:
head(ntham_rep); length(ntham_rep)
head(afr_rep); length(afr_rep)
head(nz_rep); length(nz_rep)
head(jp_rep); length(jp_rep)
head(au_rep); length(au_rep)

head(eur_al); length(eur_al)

# double check: these should all be TRUE:
table(ntham_rep %in% eur_al)
table(afr_rep %in% eur_al)
table(nz_rep %in% eur_al)
table(jp_rep %in% eur_al)
table(au_rep %in% eur_al)

ghead(eur_snp); dim(eur_snp)

propagule_result<-data.frame(iteration=1:1000, ntham_propag=NA, afr_propag=NA, nz_propag=NA, jp_propag=NA, au_propag=NA)
head(propagule_result)

regions_to_test<-data.frame(df_colname=colnames(propagule_result)[2:length(propagule_result)],allele_list=c("ntham_rep","afr_rep","nz_rep","jp_rep","au_rep"))

for (m in 1:nrow(regions_to_test)){

region.thisrun<-as.character(regions_to_test$df_colname[m])
allele.listnow<-get(as.character(regions_to_test$allele_list[m]))
head(allele.listnow); length(allele.listnow)

for (j in 1:1000){

prop_represented<-0
allele_list<-NULL
ind_list<-list()
number_propagules<-list()
i<-1

while(prop_represented<1){

number_propagules[[i]]<-paste("adding_propagule",i,sep="_")
i <- i + 1

# This line skips if the individual has already been sampled (i.e. sampling without replacement):
ind.namethisrun<-as.character(sample(eur_snp$ind,1))
if(ind.namethisrun %in% ind_list) next

ind_list<-c(ind_list,ind.namethisrun)

ind.thisrun<-unlist(eur_snp[which(eur_snp$ind %in% ind.namethisrun),3:length(eur_snp)])
head(ind.thisrun); length(ind.thisrun)
table(ind.thisrun)

alleles.thisrun<-names(ind.thisrun)[which(ind.thisrun==1)]
head(alleles.thisrun); length(alleles.thisrun)

allele_list<-unique(c(allele_list,alleles.thisrun))
head(allele_list); length(allele_list)

prop_represented<-length(allele_list)/length(allele.listnow)

} # close while

propagule_result[j, region.thisrun]<-length(number_propagules)

} # close for iterations

} # close m region

save.image("../04_workspaces/STEP06_propmin_wksp")

head(propagule_result)
tail(propagule_result)

# summarise (with replacement):
region_private$mean_prop<-apply(propagule_result[,2:length(propagule_result)],2, mean, na.rm=T)
region_private$se_prop<-apply(propagule_result[,2:length(propagule_result)],2, std.err)

region_private$n<-c(length(ntham$ind),length(afr$ind),length(nz$ind),length(jp$ind),length(au$ind))

cor.test(region_private$mean_prop, region_private$n)

region_private$unique_n<-region_private$mean_prop/region_private$n

region_private$n_rel<-region_private$unique_n*(1/min(region_private$unique_n))

# summarise (without replacement):
region_private_norepl<-region_private
region_private_norepl$mean_prop<-apply(propagule_result[,2:length(propagule_result)],2, mean, na.rm=T)
region_private_norepl$se_prop<-apply(propagule_result[,2:length(propagule_result)],2, std.err)

region_private_norepl$n<-c(length(ntham$ind),length(afr$ind),length(nz$ind),length(jp$ind),length(au$ind))

region_private_norepl$unique_n<-region_private_norepl$mean_prop/region_private_norepl$n

region_private_norepl$n_rel<-region_private_norepl$unique_n*(1/min(region_private_norepl$unique_n))

# RUN FOR 10,000

propagule_result2<-data.frame(iteration=1:10000, ntham_propag=NA, afr_propag=NA, nz_propag=NA, jp_propag=NA, au_propag=NA)
head(propagule_result2)

regions_to_test2<-data.frame(df_colname=colnames(propagule_result)[2:length(propagule_result)],allele_list=c("ntham_rep","afr_rep","nz_rep","jp_rep","au_rep"))

for (m in 1:nrow(regions_to_test2)){

region.thisrun<-as.character(regions_to_test2$df_colname[m])
allele.listnow<-get(as.character(regions_to_test2$allele_list[m]))
head(allele.listnow); length(allele.listnow)

for (j in 1:10000){

prop_represented<-0
allele_list<-NULL
number_propagules<-list()
i<-1

while(prop_represented<1){

number_propagules[[i]]<-paste("adding_propagule",i,sep="_")
i <- i + 1

ind.thisrun<-unlist(eur_snp[sample(1:nrow(eur_snp),1),3:length(eur_snp)])
head(ind.thisrun); length(ind.thisrun)
table(ind.thisrun)

alleles.thisrun<-names(ind.thisrun)[which(ind.thisrun==1)]
head(alleles.thisrun); length(alleles.thisrun)

allele_list<-unique(c(allele_list,alleles.thisrun))
head(allele_list); length(allele_list)

prop_represented<-length(allele_list)/length(allele.listnow)

} # close while

propagule_result2[j, region.thisrun]<-length(number_propagules)

} # close for iterations

} # close m region

head(propagule_result2)

# summarise:
region_private2<-region_private
region_private2$mean_prop<-apply(propagule_result2[,2:length(propagule_result2)],2, mean, na.rm=T)
region_private2$se_prop<-apply(propagule_result2[,2:length(propagule_result2)],2, std.err)

region_private2$n<-c(length(ntham$ind),length(afr$ind),length(nz$ind),length(jp$ind),length(au$ind))

region_private2$unique_n<-region_private2$mean_prop/region_private2$n

region_private2$n_rel<-region_private2$unique_n*(1/min(region_private2$unique_n))








































