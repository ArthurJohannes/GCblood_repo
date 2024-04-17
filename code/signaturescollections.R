
#################  R SCRIPT SIGNATURESCOLLECTIONS.R ###########################################

####  COMBINING AVERAGED EXPRESSION CORRELATION PROFILES FOR SINGLE GENES TO DERIVE ROBUST GENE SIGNATURE
####  USING PRECALCULATED RESULTS OBTAINED WITH R SCRIPT SELECTCOLLECTIONS.R 



###############################################################################################
####  FIRST WITH 6 GENES FOR GC SIGNATURE 1
###############################################################################################

#####################
####  read precalculated results

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagenotsevere_listall  <- readRDS ("average_correlations_notsevere.rds")
#####################

#####################
#### take 6 query genes from results

thenames <- names (averagenotsevere_listall)
querysignature1 <- c ("PER1","ZBTB16","DDIT4","TSC22D3","KLF9","IRS2")
theorder<- numeric () ; for (i in querysignature1) {thisone <- which (thenames  == i); theorder <- c (theorder, thisone)}

averagenotsevere_list <- list ()
counterb <- 0 ; for (j in theorder){counterb <- counterb +1 ;averagenotsevere_list [[counterb]] <- averagenotsevere_listall  [[j]]}
names (averagenotsevere_list) <- querysignature1
#####################

#####################
#### set use > 5 (minimal number of supporting datasets)
lengthquery1 <- length (querysignature1 )

averageprofile_dfuse <- list ()
for (i in 1:lengthquery1 ){averagedf <- averagenotsevere_list [[i]];   use <- averagedf [,3] > 5 ; averagedfuse <- averagedf [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1]; averageprofile_dfuse [[i]] <- averagedfuse }
#####################

#####################
#### combine and average 6 profiles in list averageprofile_dfuse

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:lengthquery1 )) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}

thesecolnumbers <- c (1:lengthquery1)*3 +1
merger1$meanranking <- rowMeans (merger1 [,thesecolnumbers], na.rm = TRUE)

querygenes <- names (averagenotsevere_list) [1:zoveel]
kolnames1 <- character () ; for (i in querygenes) {begin <- paste (c ("cor", "nr", "rank"), i, sep = "."); kolnames1 <- c (kolnames1, begin)}
merger1order <- merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", kolnames1,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
#####################

#####################
#### remove genes supported by too few datasets in total from all query genes

merger1order$sumofsupportingdatasets <- rowSums(merger1order [,thesecolnumbers-1], na.rm = TRUE)
hitwithmostgenequeries <- merger1order$sumofsupportingdatasets > 5 * length (querygenes)

merger2order <- merger1order [hitwithmostgenequeries,]

rownames (merger2order) <- 1:dim (merger2order) [1]
merger2order$rank.mean <- 1:dim (merger2order) [1]
#####################

#####################
#### print ranked profiles and GC signature 1 combined for 6 genes on collection of 15 datasets not severe
print ("RANKED PROFILES FOR 6 GC SIGNATURE 1 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger2order))

print (merger2order [1:50, c (1, 20)])
#####################


###############################################################################################
####  SECOND WITH 7 GENES FOR GC SIGNATURE 2
###############################################################################################

#####################
####  read precalculated results

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagesevere_listall <-  readRDS ("average_correlations_severe.rds" )
#####################

#####################
#### take 7 query genes from results

thenames <- names (averagesevere_listall )
querysignature2 <- c ("FLT3","ADORA3","CD163","OLAH","DAAM2","ADAMTS2","VSIG4")
theorder <- numeric () ; for (i in querysignature2) {thisone <- which (thenames == i); theorder <- c (theorder, thisone)}

averagesevere_list <- list ()
counterb <- 0 ; for (j in theorder){counterb <- counterb +1 ;averagesevere_list [[counterb]] <- averagesevere_listall [[j]]}
names (averagesevere_list) <- querysignature2

#####################

#####################
#### set use > 5 (minimal number of supporting datasets)

lengthquery2  <- length (querysignature2 )

thesecolnumbers <- c (1:lengthquery2)*3 +1

averageprofile_dfuse <- list ()
for (i in 1:lengthquery2){averagedf <- averagesevere_list [[i]];  use <- averagedf [,3] > 5 ; averagedfuse <- averagedf [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1]; averageprofile_dfuse [[i]] <- averagedfuse }

#####################
#### combine and average 7 profiles in list averageprofile_dfuse

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:lengthquery2 )) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}
merger1$meanranking <- rowMeans (merger1 [,thesecolnumbers], na.rm = TRUE)

querygenes <- names (averagesevere_list) [1:7]
kolnames1 <- character () ; for (i in querygenes) {begin <- paste (c ("cor", "nr", "rank"), i, sep = "."); kolnames1 <- c (kolnames1, begin)}
merger1order <- merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", kolnames1,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
merger1order$rank.mean <- 1:dim (merger1order) [1]

#####################

#####################
#### remove genes supported by too few datasets in total from all query genes

merger1order$sumofsupportingdatasets <- rowSums(merger1order [,thesecolnumbers-1], na.rm = TRUE)
hitwithmostgenequeries <- merger1order$sumofsupportingdatasets > 5 * length (querygenes)

merger2order <- merger1order [hitwithmostgenequeries,]

rownames (merger2order) <- 1:dim (merger2order) [1]
merger2order$rank.mean <- 1:dim (merger2order) [1]

#####################
#### print ranked profiles and GC signature 2 combined for 7 genes on collection of 15 datasets severe

print ("RANKED PROFILES FOR 7 GC SIGNATURE 2 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger2order))

print (merger2order [1:50, c (1, 23)])

#####################

################################################################################
#### END
################################################################################