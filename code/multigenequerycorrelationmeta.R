##  TITLE



################################################################################
#### CHOOSE MULTIGENE QUERY, AND DATASET COLLECTIONS
################################################################################


#####################
#### read all identifiers, is all genes on platforms GPL570 and GPL10558
#### here no RNAseq datasets present, otherwise add more RNAseq identifiers
#### to allow new HUGO gene names and gene aliases absent from GPL570 and GPL10558

setwd("~/GitHub/GCblood_repo/data")
identifiersGPL570 <- readRDS ("identifiers_gpl570.rds")
identifiersGPL10558 <- readRDS ("identifiers_gpl10558.rds")
identifiers <- unique (c (identifiersGPL570,identifiersGPL10558))
#####################

#####################
#### choose multigene query

#### query <- c ("TSC22D3","PER1","ZBTB16","KLF9","CXCR4","DDIT4","IRS2") 
query <- c ("ADAMTS2","CD163","VSIG4","FLT3","ADORA3","OLAH","DAAM2")
#####################

####################
## choose dataset collection

## setwd("~/GitHub/GCblood_repo/data/expressiontables_notsevere")
setwd("~/GitHub/GCblood_repo/data/expressiontables_severe")
datasetnames <- dir ()
#####################





############$$$$$$$$$$$$$$$$###################$$$$$$$$$$$$$$$$$$$



query <- intersect (query, identifiers)
querycounter <- 0
querylength <- length (query)


averageprofilesallqueries_list <- list ()
allprofilesallqueries_list <- list ()

counterdown <- length (datasetnames)*querylength

for (jj in query){

querycounter <- querycounter + 1
 
  allprofiles_list <- list ()
  ## NIEUWE REGEL 1 
  
  samplesizes <- numeric ()
  
datasetnaam <- rep (c (""), times = 1000)
gennaam  <- rep (c (""), times = 1000)
datasetnr <- rep (c (0), times = 1000)
probeaantal<- rep (c (0), times = 1000)



counterdatasetup = 1
counterprofiles = 1
for (j in datasetnames){
  
  ## counter <- counter -1 
  
  ## setwd("C:/Users/Arthur/Desktop/bloodRNA/GDStables")
  
  
  ## if (j %in% datasetnames)  {setwd("~/GitHub/GCblood_repo/data/expressiontables_notsevere"); expressiontable <- readRDS (j)}
  if (j %in% datasetnames)  {setwd("~/GitHub/GCblood_repo/data/expressiontables_severe"); expressiontable <- readRDS (j)}
  
      expressionmatrix <- t (expressiontable [, -c (1,2)])
  
  colnames (expressionmatrix) <- expressiontable [,2]
  
  probes <- which (colnames (expressionmatrix) == jj)
  
  
  if (class ( expressionmatrix [,1]) != "numeric") { expressionmatrix <- apply (expressionmatrix, 2, as.numeric)}
  
  probes <- which (colnames (expressionmatrix) == jj)
  
  
  gennaam <- jj
  aantal <- length(probes)

  
  counterprobesets = 0
  
   
  datasetnaam [counterdatasetup] <-  j
  gennaam [counterdatasetup] <- gennaam
  datasetnr [counterdatasetup] <- counterdatasetup
  probeaantal [counterdatasetup] <- aantal
  
  for (probe in probes){
  
    correlations_namednum <- cor (expressionmatrix [,probe],expressionmatrix)
    b <-  correlations_namednum [1, ]
    correlationsranked_namednum <- - sort (-b)
    df <- as.data.frame (correlationsranked_namednum)
    df$IDENTIFIER <- names (correlationsranked_namednum)
    colnames (df) [1] <-  paste (j, "nr" , probe, sep = "")
  
  ## zoveel <- counterprobesets + counterprofiles
  allprofiles_list [[counterprobesets + counterprofiles]] <- df
  
  counterprobesets <- counterprobesets +1
 
  }
  counterprofiles <- counterprofiles + counterprobesets 
  counterdown <- counterdown -1 
  counterdatasetup <- counterdatasetup +1
print ("file reads left to do:")
print (counterdown)

}



## intersect filter om beste probe per dataset te vinden

ll <- length (allprofiles_list)
matriksje <- matrix (nrow = ll, ncol = ll)
for (x in 1:ll){
  
for(y in 1:ll){
  
 matriksje [x,y]  <- length (intersect (allprofiles_list [[x]][1:200,2], allprofiles_list [[y]][1:200,2]))
}

}
matriksje[ row(matriksje) == col(matriksje) ] <- 0
kolomgemiddelde <- colMeans (matriksje)

## daarna ev filter om beste datasets te vinden (hierbij mogelijkheid van 2 veschillende hits verloren)





probeer <- rep (as.numeric (datasetnr), times = (as.numeric(probeaantal)))
beide <- cbind (probeer, kolomgemiddelde)
beide <- as.data.frame (beide)

ZOBETER <- unique (beide$probeer)

allprofilesallqueries_list [[querycounter]] <- allprofiles_list

tothier <- length (allprofiles_list)                       

vectorintallaa <- as.integer ()
for (aantaa in ZOBETER){
  aaaa <- which (beide$probeer == aantaa)
  use <-beide$kolomgemiddelde [aaaa]     >5 
  vectorintaa <- aaaa [use]
  vectorintallaa <- c (vectorintallaa, vectorintaa)
}

 





vectorintall <- as.integer ()
for (aant in ZOBETER){
aa <- which (beide$probeer == aant)
use <-beide$kolomgemiddelde [aa] == max (beide$kolomgemiddelde [aa]) 
vectorint <- aa [use] [1]
vectorintall <- c (vectorintall, vectorint)
}

## beide filters


filter <- intersect (vectorintallaa, vectorintall)
tothier2 <- filter [length (filter)]
checkfilter <- filter
if (length (filter) == 0){filter <- c (1,2)}

## weet niet waarom maar een vectorintall  [[11]] = allprofiles_list [[21,22,23]] gaat niet met fun = max, wel met mean maar geeft dan NAs
## gds 3646, daarom gebruikt vectorintallb

## allprofiles_list <- allprofiles_list
## tothier <- length (allprofiles_list)
for (j in filter){
  allprofiles_list [[j]] <- aggregate (allprofiles_list[[j]] [,1], by = list (IDENTIFIER = allprofiles_list [[j]] [,2]), FUN = max)
}
## betere selectie van groot platform met veel genes is nodig


## identifiers  <-  unique (allprofiles_list [[tothier2]] [,1])
mergedprofiles_df <- as.data.frame (identifiers) 
colnames (mergedprofiles_df) <- "IDENTIFIER"


##  omschakelen van max correlatie waarde naar rank, meancorrelation en meanrank


for (j in filter){
  
  mergedprofiles_df <- merge  (mergedprofiles_df, allprofiles_list [[j]], by = "IDENTIFIER", all.x = TRUE )
}

rlength <- nrow (mergedprofiles_df)
clength <- ncol (mergedprofiles_df)

for ( r in 1: rlength)
{mergedprofiles_df$Correlation [r] <-  rowMeans (mergedprofiles_df [r, 2:clength], na.rm = TRUE)
mergedprofiles_df$supporting_datasets [r] <- sum (!is.na(mergedprofiles_df [r, 2:clength]))

}

orderedmergedprofiles_df <- mergedprofiles_df [order (-mergedprofiles_df$Correlation), ]
averageprofiledf <- orderedmergedprofiles_df [, c (1, clength +1, clength + 2)]

rownames (averageprofiledf) <- 1: length (identifiers)

averageprofilesallqueries_list [[querycounter]] <- averageprofiledf
if (length (checkfilter) == 0) {averageprofilesallqueries_list [[querycounter]] <- "no hits above 5"}
## print (query [querycounter])
print (paste (query [querycounter], "done", sep = " "))
}
names (averageprofilesallqueries_list) <- query

############## PRINT RESULTS ##############################

for (i in 1:querylength){use <- averageprofilesallqueries_list [[i]] [,3] >5; print (averageprofilesallqueries_list [[i]] [use,][1:25,])}


## in case code above was run with aggregate set to minus to see negative correlations (tails)
##for (i in 1:querylength){use <- averageprofilesallqueries_list [[i]] [,3] >5; aa <- averageprofilesallqueries_list [[i]] [use,]; bb <- dim (aa) [1]; print (names (averageprofilesallqueries_list) [i]); cc <- bb-25 ;print (aa [cc:bb,])}



##############$$$$$$$$$$$$$$$$$$$$$$$$##################$$$$$$$$$$$$$$$$$$$$
##  hier komt code om signatures te laten zien met separate dataset collections 

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

forGCsignature1_list  <- readRDS ("average_correlations_notsevere.rds")



##   set use > 5
averageprofile_dfuse <- list ()
for (i in 1:7){vang <- forGCsignature1_list [[i]];   use <- vang [,3] > 5 ; vanguse <- vang [use,]; vanguse$ranking <- 1:dim (vanguse) [1]; averageprofile_dfuse [[i]] <- vanguse }

## merge but leave out profile for CXCR4

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3,4,6,7)) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}
merger1$meanranking <- rowMeans (merger1 [,c (4,7,10,13,16,19)], na.rm = TRUE)

## haal weg bv genes die minder dan 3 maal in 6 profiles staan, doe ook bij scan ?

querygenes <- c ("TSC22D3","PER1","ZBTB16","KLF9","DDIT4","IRS2")
allemaal <- character () ; for (i in querygenes) {begin <- paste (c ("cor", "nr", "rank"), i, sep = "."); allemaal <- c (allemaal, begin)}
merger1order <- merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", allemaal,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
merger1order$rank.mean <- 1:dim (merger1order) [1]
 
##  nu GC signature voor severe op 15, check ook in scan code dan

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

forGCsignature2_list <- readRDS ("average_correlations_severe.rds")

##   set use > 5
averageprofile_dfuse <- list ()
for (i in 1:7){vang <- forGCsignature2_list [[i]];  use <- vang [,3] > 5 ; vanguse <- vang [use,]; vanguse$ranking <- 1:dim (vanguse) [1]; averageprofile_dfuse [[i]] <- vanguse }

querygenes <- names (forGCsignature2_list)

merger1 <- merge (averageprofile_dfuse [[1]], averageprofile_dfuse [[2]], by = "IDENTIFIER", all = TRUE)
for (i in c (3:7)) {merger1 <- merge (merger1, averageprofile_dfuse [[i]], by = "IDENTIFIER", all = TRUE)}
merger1$meanranking <- rowMeans (merger1 [,c (4,7,10,13,16,19,22)], na.rm = TRUE)

## haal weg bv genes die minder dan 3 maal in 7 profiles staan, doe ook bij scan ?

allemaal <- character () ; for (i in querygenes) {begin <- paste (c ("cor", "nr", "rank"), i, sep = "."); allemaal <- c (allemaal, begin)}
merger1order <- merger1 [order (merger1$meanranking, decreasing = FALSE),]

colnames (merger1order) <- c ("gene", allemaal,"mean.rank")
rownames (merger1order) <- 1:dim (merger1order) [1]
merger1order$rank.mean <- 1:dim (merger1order) [1]














#####################################################################

setwd("~/GitHub/GCblood_repo/results/fromseparatecollections")

## averageprofilesallqueries_list  <- readRDS ("average_correlations_notsevere.rds")
averageprofilesallqueries_list  <- readRDS ("average_correlations_severe.rds" )


query <- names (averageprofilesallqueries_list)

querylengte <- length (query)


vaanglijstuse <- list ()
 ##  check nog of use groter 5 beter is dan groter 2
 
 for (i in 1:querylengte) {use <- averageprofilesallqueries_list [[i]] [,3] > 5 ; ditdf <- averageprofilesallqueries_list [[i]] [use,]; ditdf$ranking <- 1:dim (ditdf) [1]; vaanglijstuse [[i]] <- ditdf}
 ## > for (i in 1:10){print (vaanglijstuse [[i]] [1:10,])}

 mijngenes <- as.data.frame (query)
 colnames (mijngenes) [1] <- "IDENTIFIER"
 merger1 <- merge (mijngenes, vaanglijstuse [[1]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)
 for (j in 2:querylengte){ merger1 <- merge (merger1, vaanglijstuse [[j]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)}
 colnames (merger1) <- c ("gene", query)
 genesnow <- as.character (merger1$gene)
 volgorde <- numeric (); for (i in query){deze <- which (genesnow == i); volgorde <- c (volgorde, deze)}
 merger2 <- merger1 [volgorde,]
 rownames (merger2) <- as.character (merger2$gene)
 merger3 <- merger2 [,-1]
 
 ## transpose to have profiles in rows, asymmetric square matrix
 tmerger3 <- t (merger3)
 
 print (tmerger3)
 
 ###################################### dit laat voorlopig maar 
 
 ## moregenes <- c ("ZFP36L2","MCL1")
 moregenes <- c ("ARG1","MAOA")
 
 
 samengenes <- c (query, moregenes)
 
 mijngenes <- as.data.frame (samengenes)
 
 colnames (mijngenes) [1] <- "IDENTIFIER"
 merger1 <- merge (mijngenes, vaanglijstuse [[1]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)
 for (j in 2:querylengte){ merger1 <- merge (merger1, vaanglijstuse [[j]] [,c (1,4)], by = "IDENTIFIER",all.x = TRUE)}
 colnames (merger1) <- c ("gene", query)
 genesnow <- as.character (merger1$gene)
 volgorde <- numeric (); for (i in samengenes){deze <- which (genesnow == i); volgorde <- c (volgorde, deze)}
 merger2 <- merger1 [volgorde,]
 rownames (merger2) <- as.character (merger2$gene)
 merger3 <- merger2 [,-1]
 
 ##  transpose to have profiles in rows, asymmetric rectangular matrix
 tmerger3 <- t (merger3)
 
 print (tmerger3)
 
 ######################################
 
 ##  add mean ranking for gene sets or mdules i profiles 
 
 
 
 
 
 
 