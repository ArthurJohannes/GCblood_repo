
#################  R SCRIPT  SIGNATURESPROGRAMMATICALLY.R #####################################

####  COMBINING AVERAGED EXPRESSION CORRELATION PROFILES FOR SINGLE GENES TO DERIVE ROBUST GENE SIGNATURE
####  USING PRECALCULATED RESULTS OBTAINED WITH R SCRIPT SELECTPROGRAMMATICALLY.R 

###############################################################################################
####  FIRST WITH 6 GENES FOR GC SIGNATURE 1
###############################################################################################

#####################
####  read precalculated results
setwd("~/GitHub/GCblood_repo/results/signatures")

forGCsignature1_list <- readRDS ("forGCsignature1.rds")

#####################

#####################
#### set minimal number of supporting datasets required for all genes present in each averageprofiledf  : use > 5

averageprofile_dfuse <- list ()
for (i in 1:6) {
  file <-
    forGCsignature1_list [[i]]; averageprofiledf <-
      file$averageprofile_df;  use <-
        averageprofiledf [,3] > 5 ; averageprofiledfuse <-
          averageprofiledf [use,]; averageprofiledfuse$ranking <-
            1:dim (averageprofiledfuse) [1]; averageprofile_dfuse [[i]] <-
              averageprofiledfuse
}

#####################

#####################
####  combine 6 profiles in list averageprofile_dfuse for average of averages


merger1 <-
  merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3,4,5,6)) {
  merger1 <-
    merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)
}
merger1$meanranking <-
  rowMeans (merger1 [,c (4,7,10,13,16,19)], na.rm = TRUE)

querygenes <- names (forGCsignature1_list)
mergercolnames <-
  character () ; for (i in querygenes) {
    begin <-
      paste (c ("cor", "nr", "rank"), i, sep = "."); mergercolnames <-
        c (mergercolnames, begin)
  }
merger1order <-
  merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", mergercolnames,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
merger1order$rank.mean <- 1:dim (merger1order) [1]

#####################

#####################
#### remove genes supported by too few profiles in total from all query genes


merger1order$sumofsupportingprofiles <- rowSums(merger1order [,c (3,6,9,12,15,18)], na.rm = TRUE)
hitwithmostgenequeries <- merger1order$sumofsupportingprofiles > 5*length (querygenes)

merger2order <- merger1order [hitwithmostgenequeries,]

rownames (merger2order) <- 1:dim (merger2order) [1]
merger2order$rank.mean <- 1:dim (merger2order) [1]

#####################

#####################
#### print results
print ("RANKED PROFILES FOR 6 GC SIGNATURE 1 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger2order))

print (merger2order [1:50, c (1, 20)])
#####################

###############################################################################################
####  SECOND WITH 7 GENES FOR GC SIGNATURE 2
###############################################################################################

#####################
####  read precalculated results

setwd("~/GitHub/GCblood_repo/results/signatures")

forGCsignature2_list <- readRDS ("forGCsignature2.rds")
#####################

#####################
#### set minimal number of supporting datasets required for all genes present in each averageprofiledf  : use > 5

averageprofile_dfuse <- list ()
for (i in 1:7) {
  file <-
    forGCsignature2_list [[i]]; averageprofiledf <-
      file$averageprofile_df;  use <-
        averageprofiledf [,3] > 5 ; averageprofiledfuse <-
          averageprofiledf [use,]; averageprofiledfuse$ranking <-
            1:dim (averageprofiledfuse) [1]; averageprofile_dfuse [[i]] <-
              averageprofiledfuse
}

#####################

#####################
####  combine 7 profiles in list averageprofile_dfuse for average of averages

querygenes <- names (forGCsignature2_list)

merger1 <-
  merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:7)) {
  merger1 <-
    merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)
}
merger1$meanranking <-
  rowMeans (merger1 [,c (4,7,10,13,16,19,22)], na.rm = TRUE)


mergercolnames <-
  character () ; for (i in querygenes) {
    begin <-
      paste (c ("cor", "nr", "rank"), i, sep = "."); mergercolnames <-
        c (mergercolnames, begin)
  }
merger1order <-
  merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", mergercolnames,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
merger1order$rank.mean <- 1:dim (merger1order) [1]

#####################

#####################
#### remove genes supported by too few profiles in total from all query genes


merger1order$sumofsupportingprofiles <- rowSums(merger1order [,c (3,6,9,12,15,18,21)], na.rm = TRUE)
hitwithmostgenequeries <- merger1order$sumofsupportingprofiles > 5*length (querygenes)

merger2order <- merger1order [hitwithmostgenequeries,]

rownames (merger2order) <- 1:dim (merger2order) [1]
merger2order$rank.mean <- 1:dim (merger2order) [1]

#####################

#####################
#### print results

print ("RANKED PROFILES FOR 7 GC SIGNATURE 2 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger2order))

print (merger2order [1:50, c (1, 23)])

###############################################################################################
##################### END
###############################################################################################
