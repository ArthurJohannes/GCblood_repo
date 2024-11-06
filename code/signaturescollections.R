
#################  R SCRIPT SIGNATURESCOLLECTIONS.R ###########################################

####  COMBINING AVERAGED EXPRESSION CORRELATION PROFILES FOR SINGLE GENES TO DERIVE ROBUST GENE SIGNATURE
####  USING PRECALCULATED RESULTS OBTAINED WITH R SCRIPT SELECTCOLLECTIONS.R 

####  USE AT LEAST 2 GENES IN QUERY


###############################################################################################
####  FIRST WITH query GENES FOR GC SIGNATURE 1
###############################################################################################

#####################
####  read precalculated results

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagenotsevere_listall  <- readRDS ("average_correlations_notsevere.rds")
#####################

#####################
#### choose query genes for GC signature 1 (little bias, or neutrophil bias)

thenames <- names (averagenotsevere_listall)

queryGCsignature1littlebias <- c ("PER1","ZBTB16","DDIT4","TSC22D3","KLF9","IRS2")
queryGCsignature1biasneutrophil <- c ("TSC22D3","IRS2","FKBP5","ECHDC3","TPST1")


querysignature1 <- queryGCsignature1littlebias ; querychoice <- "query GC signature 1 little bias"
## querysignature1 <- queryGCsignature1biasneutrophil ; querychoice <- "query GC signature 1 bias neutrophil"

theorder<- numeric () ; for (i in querysignature1) {thisone <- which (thenames  == i); theorder <- c (theorder, thisone)}

averagenotsevere_list <- list ()
counterb <- 0 ; for (j in theorder){counterb <- counterb +1 ;averagenotsevere_list [[counterb]] <- averagenotsevere_listall  [[j]]}
names (averagenotsevere_list) <- querysignature1
#####################

#####################
#### set use > 5 (minimal number of supporting datasets for collection n = 15)
lengthquery1 <- length (querysignature1 )

averageprofile_dfuse <- list ()
for (i in 1:lengthquery1 ){averagedf <- averagenotsevere_list [[i]];   use <- averagedf [,3] > 5 ; averagedfuse <- averagedf [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1]; averageprofile_dfuse [[i]] <- averagedfuse }
#####################

#####################
#### combine and average query gene profiles in list averageprofile_dfuse

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:lengthquery1 )) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}

thesecolnumbers <- c (1:lengthquery1)*3 +1
merger1$meanranking <- rowMeans (merger1 [,thesecolnumbers], na.rm = TRUE)

querygenes <- names (averagenotsevere_list)
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
#### print ranked profiles and GC signature 1 combined for query genes on collection of 15 datasets not severe
print ("RANKED PROFILES GC SIGNATURE 1 QUERY GENES AND RANKED AVERAGED PROFILE")

print (querychoice)

print (head (merger2order))

colnrmeanrank <- which (colnames (merger2order) == "mean.rank")

print (querychoice)

print (merger2order [1:50, c (1, colnrmeanrank)])
#####################


###############################################################################################
####  SECOND WITH query GENES FOR GC SIGNATURE 2
###############################################################################################

#####################
####  read precalculated results

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagesevere_listall <-  readRDS ("average_correlations_severe.rds" )
#####################

#####################
#### choose query genes for GC signature 2

thenames <- names (averagesevere_listall )

queryGCsignature2biasmyeloid <- c ("FLT3","ADORA3","CD163","OLAH","DAAM2","ADAMTS2","VSIG4")
queryGCsignature2biasneutrophil <- c ("OLAH","IL1R2","IL18R1","FKBP5","ECHDC3")
queryGCsignature2biasmonocyte <- c ("FLT3","ADAMTS2","VSIG4","AMPH","GPER1")


querysignature2 <- queryGCsignature2biasmyeloid ; querychoice <- c ("query GC signature 2 bias myeloid cells")
## querysignature2 <- queryGCsignature2biasneutrophil ; querychoice <- c ("query GC signature 2 bias neutrophil") 
## querysignature2 <- queryGCsignature2biasmonocyte ;  querychoice <- c ("query GC signature 2 bias monocyte") 


theorder <- numeric () ; for (i in querysignature2) {thisone <- which (thenames == i); theorder <- c (theorder, thisone)}

averagesevere_list <- list ()
counterb <- 0 ; for (j in theorder){counterb <- counterb +1 ;averagesevere_list [[counterb]] <- averagesevere_listall [[j]]}
names (averagesevere_list) <- querysignature2

#####################

#####################
#### set use > 5 (minimal number of supporting datasets for collection n = 15)

lengthquery2  <- length (querysignature2 )

thesecolnumbers <- c (1:lengthquery2)*3 +1

averageprofile_dfuse <- list ()
for (i in 1:lengthquery2){averagedf <- averagesevere_list [[i]];  use <- averagedf [,3] > 5 ; averagedfuse <- averagedf [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1]; averageprofile_dfuse [[i]] <- averagedfuse }

#####################
#### combine and average query gene profiles in list averageprofile_dfuse

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:lengthquery2 )) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}
merger1$meanranking <- rowMeans (merger1 [,thesecolnumbers], na.rm = TRUE)

querygenes <- names (averagesevere_list)
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
#### print ranked profiles and GC signature 2 combined for query genes on collection of 15 datasets severe

print ("RANKED PROFILES GC SIGNATURE 2 QUERY GENES AND RANKED AVERAGED PROFILE")

print (querychoice)

print (head (merger2order))

colnrmeanrank <- which (colnames (merger2order) == "mean.rank")

print (querychoice)

print (merger2order [1:50, c (1, colnrmeanrank)])

#####################

################################################################################
#### END
################################################################################