
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
#### use 19 selected genes


thenames <- names (averagenotsevere_listall)

geneselection19 <- c ("ZBTB16","KLF9","DDIT4","TSC22D3","PER1","IRS2","CXCR4","FKBP5","CD163","VSIG4","ADAMTS2","FLT3","MAOA","AMPH","ADORA3","MACIR","OLAH","DAAM2","ARG1")

theorder<- numeric () ; for (i in geneselection19) {thisone <- which (thenames  == i); theorder <- c (theorder, thisone)}

averagenotsevere_listselect19 <- list ()
counterb <- 0 ; for (j in theorder){counterb <- counterb +1 ;averagenotsevere_listselect19 [[counterb]] <- averagenotsevere_listall  [[j]]}
names (averagenotsevere_listselect19) <- geneselection19



#####################

#####################
####  remove bad query genes without a resulting profile (gene AMPH)

listlength <- length (averagenotsevere_listselect19)

averagenotsevere_list <- list ()
goodquerynumbers <- numeric (); for (i in 1:listlength){theclass <- class (averagenotsevere_listselect19 [[i]]);if (theclass != "character"){goodquerynumbers  <- c (goodquerynumbers , i)}}
countera <- 0 ; for (i in goodquerynumbers) {countera  <- countera  +1 ; averagenotsevere_list [[countera]] <- averagenotsevere_listselect19  [[i]]}
names (averagenotsevere_list) <- names (averagenotsevere_listselect19) [goodquerynumbers]

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

print ("reciprocal gene ranking in average gene expression correlation profiles from collection WB not severe 15")

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
#### use 19 selected genes

thenames <- names (averagesevere_listall)

geneselection19 <- c ("ZBTB16","KLF9","DDIT4","TSC22D3","PER1","IRS2","CXCR4","FKBP5","CD163","VSIG4","ADAMTS2","FLT3","MAOA","AMPH","ADORA3","MACIR","OLAH","DAAM2","ARG1")

theorder<- numeric () ; for (i in geneselection19) {thisone <- which (thenames  == i); theorder <- c (theorder, thisone)}

averagesevere_listselect19 <- list ()
counterb <- 0 ; for (j in theorder){counterb <- counterb +1 ;averagesevere_listselect19 [[counterb]] <- averagesevere_listall  [[j]]}
names (averagesevere_listselect19) <- geneselection19


##################### 

#####################
####  remove bad query genes without a resulting profile (all good)
listlength <- length (averagesevere_listselect19)

averagesevere_list <- list ()
goodquerynumbers <- numeric (); for (i in 1:listlength){theclass <- class (averagesevere_listselect19[[i]]);if (theclass != "character"){goodquerynumbers  <- c (goodquerynumbers , i)}}
countera <- 0 ; for (i in goodquerynumbers) {countera  <- countera  +1 ; averagesevere_list [[countera]] <- averagesevere_listselect19 [[i]]}
names (averagesevere_list) <- names (averagesevere_listselect19) [goodquerynumbers]

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

print ("reciprocal gene ranking in average gene expression correlation profiles from collection WB severe 15")

print (reciprocal_mat)
#####################



##########################################################################
###########  END
##########################################################################





