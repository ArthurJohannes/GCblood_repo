
#################  R SCRIPT RECIPROCALCOLLECTIONS.R ################################

####################################################################################
########### RECIPROCAL RANKING OF GENES IN EXPRESSION CORRELATION PROFILES
########### USING PRECALCULATED RESULTS OBTAINED WITH R SCRIPT SELECTCOLLECTIONS.R
####################################################################################


####################################################################################
####  FIRST WITH 19 GENES INCLUDING SET GC-1 and GC-2 GENES FOR COLLECTION NOT SEVERE (N = 15)
####################################################################################

#####################
#### read precalculated results

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagenotsevere_listall  <- readRDS ("average_correlations_notsevere.rds" )

#####################

#####################
####  remove bad query genes without a resulting profile

listlength <- length (averagenotsevere_listall)

averagenotsevere_list <- list ()
goodquerynumbers <- numeric (); for (i in 1:listlength){theclass <- class (averagenotsevere_listall [[i]]);if (theclass != "character"){goodquerynumbers  <- c (goodquerynumbers , i)}}
countera <- 0 ; for (i in goodquerynumbers) {countera  <- countera  +1 ; averagenotsevere_list [[countera]] <- averagenotsevere_listall  [[i]]}
names (averagenotsevere_list) <- names (averagenotsevere_listall) [goodquerynumbers]

query <- names (averagenotsevere_list)

querylength <- length (query)

averageprofile_dfuse <- list ()

#####################

#####################
#### reciprocal gene expression correlation rankings in matrix

for (i in 1:querylength) {use <- averagenotsevere_list [[i]] [,3] > 5 ; averagedfuse <- averagenotsevere_list [[i]] [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1];averageprofile_dfuse [[i]] <- averagedfuse}


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

#####################

#####################
#### transpose to have profiles in rows, asymmetric square matrix, print matrix

reciprocal_mat <- t (merger3)

print (reciprocal_mat )
#####################

####################################################################################
####  SECOND WITH 19 GENES INCLUDING SET GC-1 and GC-2 GENES FOR COLLECTION SEVERE (N = 15)
####################################################################################

#####################
#### read precalculated results

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

averagesevere_listall  <- readRDS ("average_correlations_severe.rds" )
##################### 

#####################
####  remove bad query genes without a resulting profile
listlength <- length (averagesevere_listall)

averagesevere_list <- list ()
goodquerynumbers <- numeric (); for (i in 1:listlength){theclass <- class (averagesevere_listall [[i]]);if (theclass != "character"){goodquerynumbers  <- c (goodquerynumbers , i)}}
countera <- 0 ; for (i in goodquerynumbers) {countera  <- countera  +1 ; averagesevere_list [[countera]] <- averagesevere_listall  [[i]]}
names (averagesevere_list) <- names (averagesevere_listall) [goodquerynumbers]

query <- names (averagesevere_list)

querylength <- length (query)

averageprofile_dfuse <- list ()
#####################

#####################
#### reciprocal gene expression correlation rankings in matrix

for (i in 1:querylength) {use <- averagesevere_list [[i]] [,3] > 5 ; averagedfuse <- averagesevere_list [[i]] [use,]; averagedfuse$ranking <- 1:dim (averagedfuse) [1];averageprofile_dfuse [[i]] <- averagedfuse}


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

#####################

#####################
#### transpose to have profiles in rows, asymmetric square matrix, print matrix

reciprocal_mat <- t (merger3)

print (reciprocal_mat)
#####################



##########################################################################
###########  END
##########################################################################





