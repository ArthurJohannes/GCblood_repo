############ SHOWS GC Signatures 1 and 2  using 2 separate dataset collections
############ GC 1 Signature using multigene query 1 (n=6) on dataset collection not severe (n=15)
############ GC 2 Signature using multigene query 2 (n=7) on dataset collection severe (n=15)

#####################################################
############ PART1 COMBINING GENE PROFILES
#####################################################

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagenotsevere_list  <- readRDS ("average_correlations_notsevere.rds")



##   set use > 5
averageprofile_dfuse <- list ()
for (i in 1:7){averagedf <- averagenotsevere_list [[i]];   use <- averagedf [,3] > 5 ; averagedfuse <- averagedf [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1]; averageprofile_dfuse [[i]] <- averagedfuse }

#### using 6 query genes for GC signature 1:
#### PER1, ZBTB16, DDIT4, TSC22D3, KLF9, IRS2 

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:6)) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}
merger1$meanranking <- rowMeans (merger1 [,c (4,7,10,13,16,19)], na.rm = TRUE)
   
querygenes <- names (averagenotsevere_list) [1:6]
kolnames1 <- character () ; for (i in querygenes) {begin <- paste (c ("cor", "nr", "rank"), i, sep = "."); kolnames1 <- c (kolnames1, begin)}
merger1order <- merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", kolnames1,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]

##########
#### remove genes supported by too few datasets in total from all query genes
#### 

merger1order$sumofsupportingdatasets <- rowSums(merger1order [,c (3,6,9,12,15,18)], na.rm = TRUE)
hitwithmostgenequeries <- merger1order$sumofsupportingdatasets > 5 * length (querygenes)

merger2order <- merger1order [hitwithmostgenequeries,]

rownames (merger2order) <- 1:dim (merger2order) [1]
merger2order$rank.mean <- 1:dim (merger2order) [1]

print ("RANKED PROFILES FOR 6 GC SIGNATURE 1 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger2order))

print (merger2order [1:50, c (1, 20)])


 
##  nu GC signature voor severe op 15, check ook in scan code dan

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagesevere_list <- readRDS ("average_correlations_severe.rds")

##   set use > 5
averageprofile_dfuse <- list ()
for (i in 1:10){averagedf <- averagesevere_list [[i]];  use <- averagedf [,3] > 5 ; averagedfuse <- averagedf [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1]; averageprofile_dfuse [[i]] <- averagedfuse }


#### using 7 query genes for GC signature 2:
#### FLT3, ADORA3, CD163, OLAH, DAAM2, ADAMTS2, VSIG4


merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:7)) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}
merger1$meanranking <- rowMeans (merger1 [,c (4,7,10,13,16,19,22)], na.rm = TRUE)




querygenes <- names (averagesevere_list) [1:7]
kolnames1 <- character () ; for (i in querygenes) {begin <- paste (c ("cor", "nr", "rank"), i, sep = "."); kolnames1 <- c (kolnames1, begin)}
merger1order <- merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", kolnames1,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
merger1order$rank.mean <- 1:dim (merger1order) [1]

##########
#### remove genes supported by too few datasets in total from all query genes
#### 

merger1order$sumofsupportingdatasets <- rowSums(merger1order [,c (3,6,9,12,15,18,21)], na.rm = TRUE)
hitwithmostgenequeries <- merger1order$sumofsupportingdatasets > 5 * length (querygenes)

merger2order <- merger1order [hitwithmostgenequeries,]

rownames (merger2order) <- 1:dim (merger2order) [1]
merger2order$rank.mean <- 1:dim (merger2order) [1]

print ("RANKED PROFILES FOR 7 GC SIGNATURE 2 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger2order))

print (merger2order [1:50, c (1, 23)])






#####################################################
############ PART2 MATRIX FOR RECIPROCAL RANKING OF GENES IN PROFILES
#####################################################


########## matrix for biomarker set GC-1
########################################

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")


averagenotsevere_list  <- readRDS ("average_correlations_notsevere.rds" )


query <- names (averagenotsevere_list)

querylength <- length (query)

averageprofile_dfuse <- list ()
## vaanglijstuse <- list ()
##  check nog of use groter 5 beter is dan groter 2

for (i in 1:querylength) {use <- averagenotsevere_list [[i]] [,3] > 5 ; averagedfuse <- averagenotsevere_list [[i]] [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1];averageprofile_dfuse [[i]] <- averagedfuse}
## > for (i in 1:10){print (vaanglijstuse [[i]] [1:10,])}

querydf <- as.data.frame (query)
colnames (querydf) [1] <- "IDENTIFIER"
merger1 <- merge (querydf,averageprofile_dfuse [[1]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)
for (j in 2:querylength){ merger1 <- merge (merger1,averageprofile_dfuse [[j]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)}
colnames (merger1) <- c ("gene", query)
genesnow <- as.character (merger1$gene)
rightorder <- numeric (); for (i in query){thisone <- which (genesnow == i); rightorder <- c (rightorder, thisone)}
merger2 <- merger1 [rightorder,]
rownames (merger2) <- as.character (merger2$gene)
merger3 <- merger2 [,-1]

## transpose to have profiles in rows, asymmetric square matrix
reciprocal_mat <- t (merger3)



print (reciprocal_mat )






########## matrix for biomarker set GC-2
########################################

#####################################################################

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

## averagesevere_list  <- readRDS ("average_correlations_notsevere.rds")
averagesevere_list  <- readRDS ("average_correlations_severe.rds" )


query <- names (averagesevere_list)

querylength <- length (query)

averageprofile_dfuse <- list ()
## vaanglijstuse <- list ()
 ##  check nog of use groter 5 beter is dan groter 2
 
 for (i in 1:querylength) {use <- averagesevere_list [[i]] [,3] > 5 ; averagedfuse <- averagesevere_list [[i]] [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1];averageprofile_dfuse [[i]] <- averagedfuse}
 ## > for (i in 1:10){print (vaanglijstuse [[i]] [1:10,])}

 querydf <- as.data.frame (query)
 colnames (querydf) [1] <- "IDENTIFIER"
 merger1 <- merge (querydf,averageprofile_dfuse [[1]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)
 for (j in 2:querylength){ merger1 <- merge (merger1,averageprofile_dfuse [[j]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)}
 colnames (merger1) <- c ("gene", query)
 genesnow <- as.character (merger1$gene)
 rightorder <- numeric (); for (i in query){thisone <- which (genesnow == i); rightorder <- c (rightorder, thisone)}
 merger2 <- merger1 [rightorder,]
 rownames (merger2) <- as.character (merger2$gene)
 merger3 <- merger2 [,-1]
 
 ## transpose to have profiles in rows, asymmetric square matrix
 reciprocal_mat <- t (merger3)
 
 #### leave out CD163, not in biomarker set GC-2
 
 print (reciprocal_mat [-2,-2])
 

 
 
 #########################################
 ###########  END
 #########################################
 
 
 
 