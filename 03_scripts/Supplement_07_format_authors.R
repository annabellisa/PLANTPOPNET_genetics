# ------------------------------------ #
# ----------- SUPPLEMENT 07  --------- #
# ------------------------------------ #

# Format authors for MS

# Author: Annabel Smith

# This script streamlines formatting authors and affiliations for the MS. It doesn't do everything and it's not perfect or complete. Like, I couldn't get it to do the superscript numbers on the author names while also formatting the non-English characters. But it does assign the affiliation number to each author while maintaining special characters, which is handy. When making the numbers superscript, delete the double entries for authors with > 1 affiliation. 
# For the affiliations, I've added xx wildcards for superscripting. These can be formatted in Word using find and replace with wildcards and formats, then delete the xx after. 

library(readxl)

ath<-read_excel("../authors_final.xlsx",sheet="authors_affiliations")
ath<-as.data.frame(ath)
head(ath,3); dim(ath)

# Make sure that the table is sorted in the user-specified author order:
if(is.unsorted(ath$author_position)==T) print("need to order table by author position") else print("table is OK")

# Get list of unique affiliations in the order of the authors:
affil_list<-data.frame(affil_unique_no=1:length(unique(ath$affiliation)), affiliation=unique(ath$affiliation))
head(affil_list,3); dim(affil_list)

ath<-merge(ath, affil_list, by="affiliation", all.x=T, all.y=F)
ath<-ath[order(ath$author_position),]
ath$dumnum<-paste(ath$affil_unique_no,sep="")
head(ath); dim(ath)

auth_names<-paste(ath$name_formatted, ath$dumnum, sep="", collapse=", ")
write(auth_names, file="names.txt")

# paste affiliations together with the affiliation numbers:
affils<-paste("xx",affil_list$affil_unique_no,"xx",affil_list$affiliation,collapse=". ",sep="")

write(affils, file="affils.txt")
























