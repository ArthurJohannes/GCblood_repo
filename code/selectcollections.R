
#################  R SCRIPT SELECTCOLLECTIONS.R ################################

#### CALCULATES GENE EXPRESSION CORRELATION PROFILES AVERAGED OVER 15 DATASETS COLLECTION
#### FOR MULTIPLE GENES (PROVIDED IN A MULTIGENE QUERY) SEPARATELY

#### script suitable for both RNA seq and microarray datasets put in GEO datasets (gds) format
#### also handles multiple probes present for single gene in case of microarray data 


################################################################################
#### PART 1 : CHOOSE MULTIGENE QUERY, AND DATASET COLLECTIONS
################################################################################


#####################
#### read all gene identifiers (all gene names on platforms GPL570 and GPL10558)

setwd("~/GitHub/GCblood_repo/data")
identifiersGPL570 <- readRDS ("identifiers_gpl570.rds")
identifiersGPL10558 <- readRDS ("identifiers_gpl10558.rds")

#### to include the more recent HUGO gene names and gene aliases absent from GPL570 and GPL10558

hugoandaliasdf <- readRDS ("hugoandaliasdf.rds")
identifiersinhugoandaliasdf <- c (hugoandaliasdf$aliasname, hugoandaliasdf$hugoname)


identifiers <- unique (c (identifiersGPL570,identifiersGPL10558,identifiersinhugoandaliasdf ))
#####################

#####################
#### CHOOSE multigene query (for example either queries 1-3 with collection WB not severe,
#### or queries 4-7 with collection WB severe)

### query 1 (GC not severe little bias)
query1 <- c ("PER1","ZBTB16","DDIT4","TSC22D3","KLF9","IRS2")
### query 2 (GC not severe bias lymphoid cells)
query2 <- c ("DDIT4")
### query 3 (GC not severe bias neutrophils)
query3 <- c ("TSC22D3","IRS2","FKBP5","ECHDC3","TPST1")
### query 4 (GC severe bias lymphoid cells)
query4 <- c ("DDIT4")
### query 5 (GC severe bias myeloid cells)
query5 <- c ("FLT3","ADORA3","CD163","OLAH","DAAM2","ADAMTS2","VSIG4")
### query 6 (GC severe bias neutrophil)
query6 <- c ("OLAH","IL1R2","IL18R1","FKBP5","ECHDC3")
### query 7 (GC severe bias monocyte)
query7 <- c ("FLT3","ADAMTS2","VSIG4","AMPH","GPER1")

query19genes <- c ("ZBTB16","KLF9","DDIT4","TSC22D3","PER1","IRS2","CXCR4","FKBP5","CD163","VSIG4","ADAMTS2","FLT3","MAOA","AMPH","ADORA3","C5orf30","OLAH","DAAM2","ARG1")

##### for example query5
query <- query5
#####################

####################
#### CHOOSE dataset collections (here either collection WB not severe or collection WB severe) 
#### by setting datasetnames in lines below

setwd("~/GitHub/GCblood_repo/data/expressiontables_notsevere")
datasetnamesnotsevere <- dir ()
setwd("~/GitHub/GCblood_repo/data/expressiontables_severe")
datasetnamessevere <- dir ()

##### for example collection WB severe
## datasetnames <- datasetnamesnotsevere
datasetnames <- datasetnamessevere
#####################


################################################################################
#### PART 2 : CALCULATE AVERAGED GENE EXPRESSION CORRELATION PROFILES
################################################################################


query <- intersect (query, identifiers)
querycounter <- 0
querylength <- length (query)
allprofilesallqueries_list <- list ()
profilescore_list <- list ()
averageprofilesallqueries_list <- list ()
counterdown <- length (datasetnames)*querylength

for (querygene in query){
  
  querycounter <- querycounter + 1
  
  #### allprofiles_list will first contain gene expression correlation profiles for all probes
  #### corresponding with a single querygene, obtained for all datasets in the collection 
  
  allprofiles_list <- list ()
  samplesizes <- numeric ()
  
  ##### makes room for up to 1000 datasets in collection 
  
  namesofdatasets <- rep (c (""), times = 1000)
  genename  <- rep (c (""), times = 1000)
  datasetnr <- rep (c (0), times = 1000)
  nrofqueryprobesindatasets<- rep (c (0), times = 1000)
  
  counterdatasetup = 1
  counterprofiles = 1
  for (j in datasetnames){
    
    if (j %in% datasetnamesnotsevere)  {setwd("~/GitHub/GCblood_repo/data/expressiontables_notsevere"); expressiontable <- readRDS (j)}
    if (j %in% datasetnamessevere)  {setwd("~/GitHub/GCblood_repo/data/expressiontables_severe"); expressiontable <- readRDS (j)}
    
    #### to exchange gene aliases in table for more recent HUGO gene names
    
    genesintable <- as.character (expressiontable  [,2])
    
    genestohugo <- intersect (hugoandaliasdf$aliasname, genesintable)
    theseones <- numeric () ;for (i in genestohugo) {thisone <- which (hugoandaliasdf$aliasname == i); theseones <- c(theseones, thisone)}
    hugoandaliasdfuse <- hugoandaliasdf [theseones,]
    
    aliasestochange <- hugoandaliasdfuse$aliasname
    genesintable2 <- genesintable ;counteralias <- 0; for (i in aliasestochange) {counteralias <- counteralias +1 ; thisone <- which (genesintable == i);  genesintable2 [thisone] <- hugoandaliasdfuse$hugoname [counteralias]}
    expressiontable [,2] <- genesintable2
    
    
    
    expressionmatrix <- t (expressiontable [, -c (1,2)])
    colnames (expressionmatrix) <- expressiontable [,2]
    
    #### taking probes as rownumbers in original expression table, not as ID_REFS
    
    probes <- which (colnames (expressionmatrix) == querygene)
    
    
    if (class ( expressionmatrix [,1]) != "numeric") { expressionmatrix <- apply (expressionmatrix, 2, as.numeric)}
    
    probes <- which (colnames (expressionmatrix) == querygene)
    probenr <- length(probes)
    
    counterprobesets = 0
    
    namesofdatasets [counterdatasetup] <-  j
    genename [counterdatasetup] <- querygene
    datasetnr [counterdatasetup] <- counterdatasetup
    nrofqueryprobesindatasets [counterdatasetup] <- probenr
    
    for (probe in probes){
      
      correlations_namednum <- cor (expressionmatrix [,probe],expressionmatrix)
      b <-  correlations_namednum [1, ]
      
      #### sort to rank all gene probes in profile according to correlation value with querygene probe
      
      correlationsranked_namednum <- - sort (-b)
      df <- as.data.frame (correlationsranked_namednum)
      df$IDENTIFIER <- names (correlationsranked_namednum)
      colnames (df) [1] <-  paste (j, "nr" , probe, sep = "")
      
      
      allprofiles_list [[counterprobesets + counterprofiles]] <- df
      
      counterprobesets <- counterprobesets +1
      
    }
    counterprofiles <- counterprofiles + counterprobesets 
    counterdown <- counterdown -1 
    counterdatasetup <- counterdatasetup +1
    print ("file reads left to do:")
    print (counterdown)
    
  }
  
  
  
  #### select best correlation profiles based on consistent results of gene probes in and between datasets
  #########################################################################################################
  
  #### calculate mean intersect with other probe profiles for each separate probe profile
  #### (take gene intersects of top 200 in ranked profiles)
  
  ll <- length (allprofiles_list)
  intersectmatrix <- matrix (nrow = ll, ncol = ll)
  for (x in 1:ll){
    
    for(y in 1:ll){
      
      intersectmatrix [x,y]  <- length (intersect (allprofiles_list [[x]][1:200,2], allprofiles_list [[y]][1:200,2]))
    }
    
  }
  intersectmatrix[ row(intersectmatrix) == col(intersectmatrix) ] <- 0
  
  ####  to plot intersect matrix for visual :
  
  ##  image (intersectmatrix)
  
  kolmeans <- colMeans (intersectmatrix)
  
  #### discarding inadequate probe profiles with low mean intersect of 5 or less
  #### for smaller dataset collections use all available datasets with gene present 
  
  
  
  setnumber <- rep (as.numeric (datasetnr), times = (as.numeric(nrofqueryprobesindatasets)))
  profilescore <- cbind (setnumber, kolmeans)
  profilescore <- as.data.frame (profilescore)
  datasetnrsused <- unique (profilescore$setnumber)
  
  #### allprofilesallqueries_list will contain gene expression correlation profiles for all probes
  #### corresponding with all querygenes, obtained for all datasets in the collection 
  
  allprofilesallqueries_list [[querycounter]] <- allprofiles_list
  
  vectorintall1 <- as.integer ()
  for (xx1 in datasetnrsused){
    aa1 <- which (profilescore$setnumber == xx1)
    use <-profilescore$kolmeans [aa1]     >5 
    ##### for smaller dataset collections instead: 
    ##### use <-profilescore$kolmeans [aa1]     > 0 
    vectorint1 <- aa1 [use]
    vectorintall1 <- c (vectorintall1, vectorint1)
  }
  
  
  #### taking best probe profile for each dataset
  
  vectorintall2 <- as.integer ()
  for (xx2 in datasetnrsused){
    aa2 <- which (profilescore$setnumber == xx2)
    use <-profilescore$kolmeans [aa2] == max (profilescore$kolmeans [aa2]) 
    
    #### in case of identical maxima in one dataset, only first profile used ([1])
    
    vectorint2 <- aa2 [use] [1]
    vectorintall2 <- c (vectorintall2, vectorint2)
  }
  
  #### only keeping one best probe profile for each dataset and probe should be adequate
  #### if best probe in a dataset is inadequate, the profile (and corresponding dataset) is not used
  
  filter <- intersect (vectorintall1, vectorintall2)
  checkfilter <- filter
  
  #### keep running, in case all probes for a querygene are inadequate in all (but one) datasets, remove later
  
  checkvalues <- numeric () ; for (p in 1:ll) {checkvalue <- allprofiles_list [[p]] [1,1]; checkvalues <- c (checkvalues, checkvalue)}
  valuesone <- which (checkvalues == 1)
  
  
  if (length (filter) == 0){filter <- valuesone [1:2]}
  if (length (filter) == 1){filter <- valuesone [1:2]}
  
  
  #### for each gene in selected profiles uses the probe with maximal correlation value with the querygene
  #########################################################################################################
  
  #### CHOOSE, if especially checking for negative correlations (tails) set: FUN = min
  
  for (j in filter){
    allprofiles_list [[j]] <- aggregate (allprofiles_list[[j]] [,1], by = list (IDENTIFIER = allprofiles_list [[j]] [,2]), FUN = max)
  }
  
  
  #### use selected profiles to obtain one combined (merged) profile from all datasets with mean correlation values
  #########################################################################################################
  
  mergedprofiles_df <- as.data.frame (identifiers) 
  colnames (mergedprofiles_df) <- "IDENTIFIER"
  
  for (j in filter){
    
    mergedprofiles_df <- merge  (mergedprofiles_df, allprofiles_list [[j]], by = "IDENTIFIER", all.x = TRUE )
  }
  
  rlength <- nrow (mergedprofiles_df)
  clength <- ncol (mergedprofiles_df)
  
  for ( r in 1: rlength)
  {mergedprofiles_df$Correlation [r] <-  rowMeans (mergedprofiles_df [r, 2:clength], na.rm = TRUE)
  mergedprofiles_df$supporting_datasets [r] <- sum (!is.na(mergedprofiles_df [r, 2:clength]))
  
  }
  
  #### sort to rank all genes in combined profile according to mean correlation values with querygene
  
  orderedmergedprofiles_df <- mergedprofiles_df [order (-mergedprofiles_df$Correlation), ]
  averageprofiledf <- orderedmergedprofiles_df [, c (1, clength +1, clength + 2)]
  rownames (averageprofiledf) <- 1: length (identifiers)
  
  #### averageprofilesallqueries_list will, for each querygene, contain a single averaged gene expression correlation profile
  #### using all datasets in the collection 
  
  averageprofilesallqueries_list [[querycounter]] <- averageprofiledf
  if (length (checkfilter) == 0) {averageprofilesallqueries_list [[querycounter]] <- "no hits above threshold set at 5"}
  if (length (checkfilter) == 1) {averageprofilesallqueries_list [[querycounter]] <- "just one hit above threshold set at 5"}
  profilescore_list [[querycounter]] <- profilescore
  print (paste (query [querycounter], "done", sep = " "))
}
names (profilescore_list) <- query
names (averageprofilesallqueries_list) <- query

################################################################################
#### PART 3 : PRINT RESULTS FOR AVERAGED GENE EXPRESSION CORRELATION PROFILES
################################################################################

####  in case of query genes without a resulting profile

usegood <- logical() ; for (i in 1:querylength){useg <- class (averageprofilesallqueries_list [[i]])!= "character"; usegood <- c (usegood,useg)}
usebad <- usegood == FALSE
allquerynumbers <- 1:querylength
goodquerynumbers <- allquerynumbers [usegood]
badquerynumbers <- allquerynumbers [usebad]


#### print overview of number of sets with gene and consistent hits in sets with gene 

setnumbersunique <- numeric (); for (i in 1:querylength) {setnumberunique <- length (unique (profilescore_list [[i]] [,1])); setnumbersunique <- c (setnumbersunique, setnumberunique)}
consistenthitsall <- numeric (); for (i in 1:querylength) {if (length (intersect (i, goodquerynumbers)) == 1 ){consistenthits <- averageprofilesallqueries_list [[i]] [1,3]}else {consistenthits <- 0};consistenthitsall <- c (consistenthitsall,consistenthits ) }
overviewdf <-  as.data.frame (cbind (query, setnumbersunique, consistenthitsall)) 
colnames (overviewdf ) <- c ("gene","setswithgene","consistenthits")

print (overviewdf)

#### only shows results supported by more than 5 datasets in dataset collection (use >5)

for (i in goodquerynumbers){use <- averageprofilesallqueries_list [[i]] [,3] >5; print (averageprofilesallqueries_list [[i]] [use,][1:25,])}

#### in case code above was run with aggregate set to minus to see negative correlations (tails):
##   for (i in goodquerynumbers){use <- averageprofilesallqueries_list [[i]] [,3] >5; aa <- averageprofilesallqueries_list [[i]] [use,]; bb <- dim (aa) [1]; print (names (averageprofilesallqueries_list) [i]); cc <- bb-25 ;print (aa [bb:cc,])}

for (i in badquerynumbers) {print (query [i]); print (head (averageprofilesallqueries_list [[i]]))}

################################################################################
#### END
################################################################################



