
# ------------------------------------ #
# ------------- STEP 01  ------------- #
# ------------------------------------ #

### process SNP genotypes
### Author: Annabel Smith

# Load workspace:
load("../04_workspaces/STEP01_proc_wksp")

# load functions:
invisible(lapply(paste("../02_analysis_libraries/",dir("../02_analysis_libraries"),sep=""),function(x) source(x)))

#########################################
####     IMPORT & FORMAT DATA:    	 ####
#########################################
{

data.dir<-"../GENOTYPE_processing/DATA_FILES/FULL_DATA_SET"
dir(data.dir)

### ~~~ SINGLE ROW SNP data ~~~ ###

# The first six rows contain metadata:
gt_onerow<-read.csv(paste(data.dir,"Report_DPlan18-3451_SNP_singlerow.csv",sep="/"),header=T,skip=6)
head(gt_onerow[,1:20])

### ---- LOCUS INFO ---- ###

# First 17 columns have locus info:
loc_info_all<-colnames(gt_onerow[,1:17])

linf<-gt_onerow[,loc_info_all]
linf<-tidy.df(linf)
linf<-data.frame(locus=paste("L",1:nrow(linf),sep=""),linf)
linf$locus<-as.character(linf$locus)
head(linf)

### ---- END LOCUS INFO ---- ###

# Indices for SNPs cols, containing individual genotypes:
snp_cols<-colnames(gt_onerow)[-which(colnames(gt_onerow) %in% loc_info_all)]
head(snp_cols)
length(snp_cols)

# Make locus columns characters:
gt_onerow[,snp_cols]<-apply(gt_onerow[,snp_cols],2,as.character)

# Replace - with NA:
gt_onerow[,snp_cols]<-as.data.frame(apply(gt_onerow[,snp_cols],2,function(x)gsub("-","NA",x)))

# Make locus columns numeric:
# (gives warning: NAs introduced by coercion)
gt_onerow[,snp_cols]<-apply(gt_onerow[,snp_cols],2,as.numeric)
head(gt_onerow[,1:20])

# SNP data with individuals on rows:
snp_onerow<-data.frame(t(gt_onerow[,(which(colnames(gt_onerow)=="RepAvg")+1):length(gt_onerow)]))
colnames(snp_onerow)<-linf$locus
ind.names<-rownames(snp_onerow)

# Add individuals names back in:
snp_onerow<-data.frame(ind=ind.names,snp_onerow)
ghead(snp_onerow)

# Add site:
snp_onerow$site<-substr(snp_onerow$ind,1,unlist(lapply(gregexpr("[A-Z]", snp_onerow$ind), function(x) max(x))))
snp_onerow<-tidy.df(snp_onerow)

# Reorder:
snp_onerow<-snp_onerow[,c(which(colnames(snp_onerow) %in% c("site","ind"))[c(2,1)],grep("L",colnames(snp_onerow)))]
snp_onerow<-snp_onerow[order(snp_onerow$site,snp_onerow$ind),]
snp_onerow<-tidy.df(snp_onerow)
snp_onerow$site<-as.factor(snp_onerow$site)
ghead(snp_onerow)

### ~~~ MANUALLY FIX SITE COLS ~~~ ###
snp_onerow[,1:2]

# MAKE SITE CHARACTER TO ALLOW CHANGES:
snp_onerow$site<-as.character(snp_onerow$site)

# AL >> AL1
snp_onerow$site[snp_onerow$site=="AL"]<-"AL1"

# The original LK site was included in the new data. The original is "LK" while the new ones all had the underscore. 
# Remove LK rows without underscore:
snp_onerow[grep("LK",snp_onerow$site),1:10]

lk_rows<-which(snp_onerow$site=="LK")
old_lk_rows<-lk_rows[which(lk_rows %in% lk_rows[grep("_",snp_onerow[lk_rows,]$ind)]==F)]
snp_onerow<-snp_onerow[-old_lk_rows,]

# Change site to LK1 to match site data:
snp_onerow$site[snp_onerow$site=="LK"]<-"LK1"
snp_onerow<-tidy.df(snp_onerow)
ghead(snp_onerow)

# Change VIR to VA to match site data:
snp_onerow[grep("VIR",snp_onerow$site),1:10]
snp_onerow$site[snp_onerow$site=="VIR"]<-"VA"
snp_onerow<-tidy.df(snp_onerow)
ghead(snp_onerow)

# The outgroups are labelled "OGC" or "OGM" and for "OGM" there are more than one group. All OGs include an underscore.
og_rows<-grep("OG",snp_onerow$site)
snp_onerow[og_rows,1:10]
snp_onerow$site[og_rows]<-substr(as.character(snp_onerow$ind[og_rows]),1,(unlist(gregexpr("_",as.character(snp_onerow$ind[og_rows])))-1))

# All of the swiss sites include underscores. Their full site name occurs before the underscore.
sw_rows<-grep("SW",snp_onerow$site)
snp_onerow[sw_rows,1:10]
snp_onerow$site[sw_rows]<-substr(as.character(snp_onerow$ind[sw_rows]),1,(unlist(gregexpr("_",as.character(snp_onerow$ind[sw_rows])))-1))

# All greek sites include underscores. Their full site name occurs before the underscore.
gr_rows<-grep("GR",snp_onerow$site)
snp_onerow[gr_rows,1:10]
snp_onerow$site[gr_rows]<-substr(as.character(snp_onerow$ind[gr_rows]),1,(unlist(gregexpr("_",as.character(snp_onerow$ind[gr_rows])))-1))

# RE-FACTORISE AND TIDY:
snp_onerow$site<-as.factor(snp_onerow$site)
snp_onerow<-tidy.df(snp_onerow)
ghead(snp_onerow)
levels(snp_onerow$site)

### ~~~ END FIX SITE COLS ~~~ ###

# Import site data:
sdat<-read.table(paste("../01_data/site_data.txt",sep=""),header=T)
sdat<-sdat[sdat$genetics=="Y",]
sdat<-sdat[-which(sdat$site_code %in% c("RCH")),]
sdat<-tidy.df(sdat)
head(sdat)

# The only site not in site data should be RCH and FS (firescape):
levels(snp_onerow$site)[which(levels(snp_onerow$site) %in% sdat$site_code==F)]

snp_onerow[,1:2]

} # close import

ghead(snp_onerow)
dim(snp_onerow)
head(sdat)
dim(sdat)
levels(snp_onerow$site)

# The result from IMPORT & FORMAT DATA includes ALL 75376 loci; all sites; RCH is included because it has dups which need to be analysed and excluded separately. FS is also included. 

# 641 individuals including dups, two columns for site and individual. 

#########################################
#  Check & remove duplicate genotypes   #
#  		Save & remove RCH & FS  	    #
#########################################
{
ghead(snp_onerow)
snp_onerow[,1:2]

# Remove dups and save file with all duplicate samples for checking. All dups indicated by ".":
dup_samples<-unique(substr(snp_onerow$ind[grep("\\.",snp_onerow$ind)],1,unlist(gregexpr("\\.",snp_onerow$ind[grep("\\.",snp_onerow$ind)]))-1))
ind_all<-snp_onerow$ind
ind_all[grep("\\.",snp_onerow$ind)]<-substr(snp_onerow$ind[grep("\\.",snp_onerow$ind)],1,unlist(gregexpr("\\.",snp_onerow$ind[grep("\\.",snp_onerow$ind)]))-1)
# check
# data.frame(ind_all,snp_onerow$ind)

# Create df for all dups:
dup_data<-snp_onerow[ind_all %in% dup_samples,]
dup_data<-tidy.df(dup_data)
ghead(dup_data)

# Remove dups from main data:
snp_onerow<-snp_onerow[-grep("\\.",snp_onerow$ind),]
snp_onerow<-tidy.df(snp_onerow)
ghead(snp_onerow)

# Remove and save RCH for later analysis
rch_onerow<-snp_onerow[which(snp_onerow$site=="RCH"),]
rch_onerow<-tidy.df(rch_onerow)
ghead(rch_onerow)

snp_onerow<-snp_onerow[-which(snp_onerow$site=="RCH"),]
snp_onerow<-tidy.df(snp_onerow)
ghead(snp_onerow)

# Remove and save FS for later analysis
fs_onerow<-snp_onerow[which(snp_onerow$site=="FS"),]
fs_onerow<-tidy.df(fs_onerow)
ghead(fs_onerow)

snp_onerow<-snp_onerow[-which(snp_onerow$site=="FS"),]
snp_onerow<-tidy.df(snp_onerow)
ghead(snp_onerow)

dloci<-colnames(dup_data)[grep("L",colnames(dup_data))]
ghead(dup_data)

# What percentage of loci do not match in the duplicate data? 

dup_inds<-unique(levels(dup_data$site))
dup_res<-data.frame(ind=dup_inds,mistyped=NA,correct_typed=NA)

for (i in 1:length(dup_inds)){

ind.thisrun<-dup_inds[i]

data.thisrun<-dup_data[grep(ind.thisrun,dup_data$ind),-as.numeric(unlist(apply(dup_data[grep(ind.thisrun,dup_data$ind),3:length(dup_data)],1,function(x)which(is.na(x))))+2)]
data.thisrun[,c(1,2,sample(3:length(data.thisrun),size=10))]
ghead(data.thisrun)

t.thisrun<-table(apply(data.thisrun[,3:length(data.thisrun)],2,function(x)length(unique(x))))

dup_res[i,2]<-t.thisrun[2]
dup_res[i,3]<-t.thisrun[1]

}

dup_res$prop_mistyped<-dup_res$mistyped/(dup_res$mistyped+dup_res$correct_typed)
dup_res
head(linf) 

} # close dups

#########################################
####     RESULT = FULL DATA SET:   	 ####
#########################################
{
# The full data set includes ALL 75376 loci; all sites except RCH & FS. No dups. 

# Should be 513 individuals, 75376 loci, plus two columns for site and individual
ghead(snp_onerow)
dim(snp_onerow)

# RCH data includes all 75376 loci, 23 genotypes; no dups. Use for determining the clonality:
ghead(rch_onerow)
dim(rch_onerow)

# FS data includes all 75376 loci, 92 genotypes; no dups. 
ghead(fs_onerow)
dim(fs_onerow)

# dup_data includes 22 genotypes for seven indiviudals and was used in the last step to assess reproducibility (although these metrics were also give by dartseq and are used in the filters in the next stage)
ghead(dup_data)
dim(dup_data)

# onerow == alleles coded using dartseq notation:
# 0 = Reference allele homozygote
# 1 = SNP allele homozygote
# 2 = heterozygote

# It is wrong to code presence/absences as alleles - they indicate the opposite of the true genotype. E.g. 0,1 indicates that only the SNP allele is present, i.e. a SNP homozygote; 1,1 indicates both alleles are present, i.e. a heterozygote. If these were coded as alleles 0,1 would look like a heterozygot and 1,1 a homozygote but the opposite is true. As such, this is critical.

dim(snp_onerow)
ghead(snp_onerow)

# Site data includes all 65 sites (i.e. all except RCH & FS):
head(sdat)

# All locus info:
head(linf)

} # close full data

#########################################
####  IDENTIFY DUPLICATE SEQUENCES   ####
####   	AND FORMAT LOCUS DATA:  	 ####
#########################################
{

# Things to consider:

# There are duplicated sequences; these can be removed, even if they're different SNPs, to reduce linkage
# Even after removing duplicated sequences, there are duplicated loci because of single base pair differences (e.g. the same AlleleID but one bp difference in the sequence because of genotyping error)
# Not all duplicated sequences have the same AlleleID, so duplicated TrimmedSequences AND duplicated AlleleIDs will need to be removed
# However, even after these steps are taken, we still have duplicates. Sometimes a SNP was scored in one direction and again in the other. E.g. it was scored as a G>A SNP on one row and a A>G SNP on another row. These will, by definition, have different sequences, even if there are no genotyping errors. They can also have different AlleleIDs. So use blastn to get a similarity score for each locus. 
# If we simply remove the second record of a duplicated sequence or AlleleID, we will lose information, because there are many cases where the first has a very low call rate and the second has a much higher call rate. So, before removing duplicates, the loci will need to be ordered by call rate, or some other quality metric, to ensure the best ones are kept. 

head(linf[,1:7],2)

# some programs (e.g. bowtie2) read the > character as a separator, so update the allele ID so that the sequences can be identified in the output:
seq2<-linf[,c("AlleleID","TrimmedSequence")]
seq2$AlleleID<-as.character(seq2$AlleleID)
seq2$TrimmedSequence<-as.character(seq2$TrimmedSequence)
seq2$ai2<-seq2$AlleleID

# get first part of the seq ref for the ID:
seq2$ai3<-substr(seq2$ai2,1,unlist(gregexpr("F",seq2$ai2))+1)
seq2$ai3<-substr(seq2$ai3,1,nchar(seq2$ai3)-1)
seq2$ai3<-substr(seq2$ai3,1,nchar(seq2$ai3)-2)
head(seq2)

# get SNP position
seq2$pos<-lapply(seq2$ai2,function(x)substr(x,unlist(gregexpr("-",x))[[1]]+1,unlist(gregexpr(":",x))[[1]]-1))

# get SNP type:
seq2$snp1<-substr(seq2$ai2,(nchar(seq2$ai2)-2),(nchar(seq2$ai2)-2))
seq2$snp2<-substr(seq2$ai2,(nchar(seq2$ai2)),(nchar(seq2$ai2)))
seq2$snp<-paste(seq2$snp1,seq2$snp2,sep="_")

# And make new allele ID:
seq2$ai4<-paste(seq2$ai3,seq2$pos,seq2$snp,sep="_")
head(seq2)
check.rows(seq2)

# Some allele IDs are duplicated (see notes in first step) which will be fixed. However, some software (e.g. blast) won't run with duplicated allele IDs, so we need to make these unique:
seq2$ai4[which(duplicated(seq2$ai4))]
seq2$ai5<-paste(seq2$ai4,"_a",sep="")

# If none go over two, it's OK to use the which(duplicated()) index:
range(as.numeric(table(seq2$ai4)))
seq2$ai5[which(duplicated(seq2$ai5))]<-paste(seq2$ai4[which(duplicated(seq2$ai4))],"_b",sep="")

# This should be zero:
length(which(duplicated(seq2$ai5)))

# Check:
seq2[which(seq2$ai4==seq2$ai4[sample(which(duplicated(seq2$ai4)),1)]),]
seq2[which(seq2$ai4==seq2$ai4[sample(which(!duplicated(seq2$ai4)),1)]),]

# Write fasta file with all sequences:

for (i in 1:nrow(seq2)){

data.thisrun<-seq2[i,]

write(paste(">",data.thisrun$ai5,sep=""),file="DPlan18_seqs.fa",append=T)
write(data.thisrun$TrimmedSequence,file="DPlan18_seqs.fa",append=T)

} # close for

# Import blast results:
blast_dir<-"../GENOTYPE_processing/ANALYSIS_RESULTS/BLAST"
dir(blast_dir)

# Select lines that bound data:
start_line<-grep("hits",readLines(paste(blast_dir,"DPlan18_result.out",sep="/")))+1
end_line<-grep("BLASTN 2.7.1+",readLines(paste(blast_dir,"DPlan18_result.out",sep="/")))-1
end_line<-c(end_line[2:length(end_line)],(length(readLines(paste(blast_dir,"DPlan18_result.out",sep="/")))-1))
lines_to_read<-unlist(mapply(seq,from=start_line,to=end_line))

# Read data lines:
blast_res<-strsplit(readLines(paste(blast_dir,"DPlan18_result.out",sep="/"))[lines_to_read],split="\t")
blast_res<-data.frame(do.call(rbind,blast_res))
colnames(blast_res)<-c("query_acc_ver", "subject_acc_ver", "perc_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue"," bit_score")
# gives NA warning
blast_res<-data.frame(blast_res[,1:2],apply(blast_res[,3:ncol(blast_res)],2,function(x) as.numeric(as.character(x))))

# Remove rows where the match is with itself:
blast_res<-blast_res[,1:3]
blast_res$same<-ifelse(as.character(blast_res$query_acc_ver)==as.character(blast_res$subject_acc_ver),1,0)
blast_res<-blast_res[which(blast_res$same==0),]

# linf and seq2 have not been ordered, so if their AlleleIDs still match, the new ai5 can be added directly to linf:
dim(linf)[1]==dim(seq2)[1]
table(linf$AlleleID==seq2$AlleleID)

linf<-data.frame(locus=linf$locus,ai5=seq2$ai5,linf[,which(colnames(linf)=="AlleleID"):ncol(linf)])
head(linf,2)

length(which(duplicated(linf$AlleleID)))
length(which(duplicated(linf$ai5)))


# Order by call rate, so the best of the duplicates will be kept:
linf2<-linf[order(-linf$CallRate),]
linf2<-tidy.df(linf2)

# Identify direct duplicate sequences:
dup_seq<-as.character(linf2$locus[which(duplicated(linf2$TrimmedSequence))])
head(dup_seq)

# Identify duplicate AlleleID:
# Some remain from single bp differences
dup_aid<-as.character(linf2$locus[which(duplicated(linf2$AlleleID))])
head(dup_aid)

# Update linf2, removing the straight-forward dups:
strfwd<-unique(c(dup_seq,dup_aid))
linf2<-linf2[-which(linf2$locus %in% strfwd),]
linf2<-tidy.df(linf2)

# (3) Remove significant matches from blast:
blast_res$query_acc_ver<-as.character(blast_res$query_acc_ver)
blast_res$subject_acc_ver<-as.character(blast_res$subject_acc_ver)
all.matchloci<-unique(c(as.character(blast_res$query_acc_ver),as.character(blast_res$subject_acc_ver)))

# Many of the blast results have already been removed through dup sequences and AlleleIDs, so can simplify the blast results:
all.matchloci<-all.matchloci[which(all.matchloci %in% linf2$ai5)]
blast_res<-blast_res[unique(c(which(blast_res$query_acc_ver %in% all.matchloci),which(blast_res$subject_acc_ver %in% all.matchloci))),]
blast_res<-tidy.df(blast_res)

removed.loci<-list()

# Very slow, several minutes

for (i in 1:nrow(linf2)){

ai.thisrun<-as.character(linf2$ai5[i])

# If it's unique, it won't appear in the blast results; skip to next
if(ai.thisrun %in% all.matchloci==F) next

if(length(which(blast_res$query_acc_ver %in% ai.thisrun))>0) lines_i<-which(blast_res$query_acc_ver %in% ai.thisrun)
if(length(which(blast_res$subject_acc_ver %in% ai.thisrun))>0) lines_j<-which(blast_res$subject_acc_ver %in% ai.thisrun)

blast.thisrun<-blast_res[c(lines_i,lines_j),]
matches<-unique(c(blast.thisrun$query_acc_ver, blast.thisrun$subject_acc_ver))

# The data containing matches:
linf.thisrun<-linf2[which(linf2$ai5 %in% matches),]

# Sometimes the match(es) will already have been removed, so skip if this is the case:
if(nrow(linf.thisrun)==1) next

if(nrow(linf.thisrun)>1) print(paste("at i = ",i," there are ",nrow(linf.thisrun)," lines; performing removal",sep=""))

# This is ordered by call rate, so remove all of the loci that appear below the first
l.toremove<-as.character(linf.thisrun$ai5[2:nrow(linf.thisrun)])

# Which of the loci to remove have not already been added?
if(length(which(l.toremove %in% unlist(removed.loci)==F))>0) l.toremove<-l.toremove[which(l.toremove %in% unlist(removed.loci)==F)] else l.toremove<-NULL

if(length(l.toremove)>0) removed.loci[[i]]<-l.toremove

if(length(l.toremove)>0) print (paste("removing the following loci: ", l.toremove, sep=""))

} # close for 

save.image("../04_workspaces/STEP01_proc_wksp")

loc_res<-unlist(removed.loci)
head(loc_res,10); length(loc_res)
dup_blast<-as.character(linf2$locus[linf2$ai5 %in% loc_res])

# Add all dups to the main locus info data frame:
head(dup_seq)
head(dup_aid)
head(dup_blast)

all_dups<-unique(c(dup_seq,dup_aid,dup_blast))
linf$duplicate<-ifelse(linf$locus %in% all_dups,1,0)
head(linf,2)
check.rows(linf)


} # close sequence and locus format







