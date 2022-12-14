####  selecting transcriptomic datasets with strongly correlated query gene and module expression
####  obtaining gene expression correlation profiles, for one query gene, and average for module.


################################################################################
#### CHOOSE GENE QUERY , DATASETS, AND TESTGENES
################################################################################

#### read all identifiers, is all genes on platforms GPL570 and GPL10558
#### here no RNAseq datasets present, otherwise add more RNAseq identifiers
#### to allow new HUGO gene names and gene aliases absent from GPL570 and GPL10558


setwd("~/GitHub/GCblood_repo/data")

identifiersGPL570 <- readRDS ("identifiers_gpl570.rds")

identifiersGPL10558 <- readRDS ("identifiers_gpl10558.rds")

identifiers <- unique (c (identifiersGPL570,identifiersGPL10558))
#####################

#####################
#### choose single query gene

## query <- "TSC22D3"
query <- "ADAMTS2"

query <- intersect (query, identifiers)

if (length (query) == 0) {
  print ("gene absent, try gene alias")
}
#####################

#####################
#### choose testgenes (gene module)

## testgenes <- c ("TSC22D3","PER1","ZBTB16","KLF9","CXCR4","DDIT4","IRS2")
testgenes <-
  c ("ADAMTS2","CD163","VSIG4","FLT3","ADORA3","OLAH","DAAM2")
#####################

#####################
#### read file names of dataset expression tables to use
setwd("~/GitHub/GCblood_repo/data/expressiontables_notsevere")

datasetnamesnotsevere <- dir ()

setwd("~/GitHub/GCblood_repo/data/expressiontables_severe")

datasetnamessevere <- dir ()

datasetnames <- c (datasetnamesnotsevere,datasetnamessevere)


##########################################################################################
####  collect gene expression profiles with query for each probe and in each dataset in allprofiles_list
##########################################################################################


allprofiles_list <- list ()


#################
## just to keep some infos, room for up to 1000 profiles
namesofdatasets <- rep (c (""), times = 1000)
datasetnr <- rep (c (0), times = 1000)
nrofqueryprobesindataset <- rep (c (0), times = 1000)
#################

counterdatasetsdown <- length (datasetnames)
counterdatasetsup = 1
counterprofiles = 1

for (j in datasetnames) {
  if (j %in% datasetnamesnotsevere)  {
    setwd("~/GitHub/GCblood_repo/data/expressiontables_notsevere"); expression_df <-
      readRDS (j)
  }
  
  if (j %in% datasetnamessevere)  {
    setwd("~/GitHub/GCblood_repo/data/expressiontables_severe"); expression_df <-
      readRDS (j)
  }
  
  
  
  
  expression_mat <- t (expression_df [,-c (1,2)])
  
  colnames (expression_mat) <- expression_df [,2]
  
  
  probes <- which (colnames (expression_mat) == query)
  
  
  if (class (expression_mat [,1]) != "numeric") {
    expression_mat <- apply (expression_mat, 2, as.numeric)
  }
  
  probes <- which (colnames (expression_mat) == query)
  
  
  ######################
  ## just to keep some infos
  nrofprobesets <- length(probes)
  namesofdatasets [counterdatasetsup] <-  j
  datasetnr [counterdatasetsup] <- counterdatasetsup
  nrofqueryprobesindataset [counterdatasetsup] <- nrofprobesets
  ######################
  
  counterprobesets = 0
  
  for (probe in probes) {
    correlations_mat <- cor (expression_mat [,probe],expression_mat)
    correlations_namednum <-  correlations_mat [1,]
    correlationsranked_namednum <- -sort (-correlations_namednum)
    df <- as.data.frame (correlationsranked_namednum)
    df$IDENTIFIER <- names (correlationsranked_namednum)
    colnames (df) [1] <-  paste (j, "nr" , probe, sep = "")
    
    
    allprofiles_list [[counterprobesets + counterprofiles]] <- df
    
    counterprobesets <- counterprobesets + 1
    
  }
  counterprofiles <- counterprofiles + counterprobesets
  counterdatasetsdown <- counterdatasetsdown - 1
  counterdatasetsup <- counterdatasetsup + 1
  print ("files left to do:")
  print (counterdatasetsdown)
  
}

###############################################################################
#### Ranking gene expression profiles based on nr of testgenes hits
###############################################################################


#######################
## set number of top ranking genes in profile used for determining intersect with testgenes (now 25)
forintersect <- 25
#######################

allprofiles_listlength <- length (allprofiles_list)
profilename <-
  character (); profilenr <-
  numeric (); intersectsize <-
  numeric () ;intersectgene_list <-
  list();for (i  in 1:allprofiles_listlength) {
    aa <-
      allprofiles_list [[i]]; bb <-
        as.character (aa [1:forintersect,2]); cc <-
          intersect (testgenes,bb);profilename <-
            c (profilename, colnames (aa) [1]); profilenr <-
              c (profilenr, i); intersectsize <-
                c (intersectsize,length (cc)); intersectgene_list [[i]] <- cc
  }
df <- as.data.frame (cbind (profilename,profilenr,intersectsize))
df$profilename <- as.character(df$profilename)
df$profilenr <- as.numeric (as.character(df$profilenr))
df$intersectsize <- as.numeric (as.character(df$intersectsize))
dforder1 <- df [order (df$intersectsize,decreasing = TRUE),]

#### also keeping info for testgenes present in intersects, listing in right order

rightorder <- as.numeric (rownames (dforder1))
lengthrightorder <- length (rightorder)
intersectgeneordered_list <- list ()
for (i in 1:lengthrightorder) {
  intersectgeneordered_list [[i]] <-
    intersectgene_list [[rightorder [i]]]
}


rownames (dforder1) <- 1:dim (dforder1) [1]



###############################################################################
#### selecting top 20 profiles, aggregating genes for maximal correlation value
#### and saving in selectedprofiles_list
###############################################################################

bestprofilenrs <- dforder1$profilenr [1:20]

selectedprofiles_list <- list ()

counternew <- 0
for (j in bestprofilenrs) {
  counternew <-
    counternew + 1 ;selectedprofiles_list [[counternew]] <-
      aggregate (
        allprofiles_list[[j]] [,1], by = list (IDENTIFIER = allprofiles_list [[j]] [,2]), FUN = max
      )
}
###############################################################################
#### averaging selected gene expression correlation profiles"
###############################################################################


averageprofile_df <- list ()

mergedprofiles_df <- as.data.frame (identifiers)
colnames (mergedprofiles_df) <- "IDENTIFIER"

lengthbestprofilenrs <- length (bestprofilenrs)



for (j in 1:lengthbestprofilenrs) {
  mergedprofiles_df <-
    merge  (mergedprofiles_df, selectedprofiles_list [[j]], by = "IDENTIFIER", all.x = TRUE)
}

rlength <- nrow (mergedprofiles_df)
clength <- ncol (mergedprofiles_df)


for (r in 1:rlength)
{
  mergedprofiles_df$Correlation [r] <-
    rowMeans (mergedprofiles_df [r, 2:clength], na.rm = TRUE)
  mergedprofiles_df$supporting_profiles [r] <-
    sum (!is.na(mergedprofiles_df [r, 2:clength]))
  
}

orderedmergedprofiles_df <-
  mergedprofiles_df [order (-mergedprofiles_df$Correlation),]
averageprofile_df <-
  orderedmergedprofiles_df [, c (1, clength + 1, clength + 2)]

rownames (averageprofile_df) <- 1:length (identifiers)

## saving final results and settings in  Results list,
## but not large allprofiles_list and selectedprofiles_list


Results <- list ()

Results$alldatasetnames <- datasetnames
Results$query <- query
Results$testgenes <- testgenes
Results$topgenesnr <- length (bb)
Results$intersectgenes <- intersectgeneordered_list
Results$rankedprofiles <- dforder1
Results$selectedprofilenrs <- bestprofilenrs
Results$averageprofile_df <- averageprofile_df

names (Results$intersectgenes) <-
  as.character (Results$rankedprofiles$profilename)



#### check distribution of intersect sizes, optional, go back to inputs
#### "forintersect" (now 25: line 144) and "bestprofilenrs" (now 1:20: line 187) to optimize
#### and rerun

print ("RANKED CORRELATION PROFILES FOR QUERY GENE ADAMTS2 FROM 30 DATASETS")

print (names (Results))
print (Results$rankedprofiles)

################################
## set minimally required number of supporting profiles for gene in average profile (now 6)

nrofsupportingprofiles <- 6
################################

use <-  Results$averageprofile_df [,3] > nrofsupportingprofiles - 1

print ("AVERAGED RANKED CORRELATION PROFILE FOR QUERY GENE ADAMTS2 FROM 20 PROFILES")

print (Results$averageprofile_df [use,] [1:50,])



###################################################################################
#### USING PRECALCULATED averaged profiles on same 30 datasets as used for ADAMTS2
###################################################################################


##  FIRST GC-1

setwd("~/GitHub/GCblood_repo/results/signatures")

forGCsignature1_list <- readRDS ("forGCsignature1.rds")

##   set minimal number of supporting datasets required for all genes present in each averageprofiledf  : use > 5

averageprofile_dfuse <- list ()
for (i in 1:7) {
  file <-
    forGCsignature1_list [[i]]; averageprofiledf <-
      file$averageprofile_df;  use <-
        averageprofiledf [,3] > 5 ; averageprofiledfuse <-
          averageprofiledf [use,]; averageprofiledfuse$ranking <-
            1:dim (averageprofiledfuse) [1]; averageprofile_dfuse [[i]] <-
              averageprofiledfuse
}

##  combine 6 averageprofiledf for average of averages
##  but leave out profile for CXCR4, according to signature GC-1 in manuscript

merger1 <-
  merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3,4,6,7)) {
  merger1 <-
    merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)
}
merger1$meanranking <-
  rowMeans (merger1 [,c (4,7,10,13,16,19)], na.rm = TRUE)

querygenes <- c ("TSC22D3","PER1","ZBTB16","KLF9","DDIT4","IRS2")
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


print ("RANKED PROFILES FOR 6 GC-1 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger1order))

print (merger1order [1:50, c (1, 20)])


##  SECOND GC2

setwd("~/GitHub/GCblood_repo/results/signatures")

forGCsignature2_list <- readRDS ("forGCsignature2.rds")

## set minimal number of supporting datasets required for all genes present in each averageprofiledf  : use > 5

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

##  combine 7 averageprofiledf for average of averages


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

print ("RANKED PROFILES FOR 7 GC-2 QUERY GENES AND RANKED AVERAGED PROFILE")

print (head (merger1order))

print (merger1order [1:50, c (1, 23)])

######################################################################################
##################### END
######################################################################################
