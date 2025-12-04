####### PLOT GENE (SETS) EXPRESSION VALUE DISTRIBUTION FOR MICROARRAY GSE68310 SEASONAL RESP. VIR. INFECTION ##

library (gridExtra)
library (ggplot2)


##  read gene expression table, metadata, and gene modules to use

setwd("~/GitHub/GCblood_repo/plotting/plotdata")

## read table gse68310 pre-subsetted for 22 genes  (file size limit)
gse68310 <- readRDS ("Tablegse68310sub22genes")
gse68310meta <- readRDS ("gse68310meta.rds")
modules <- readRDS ("modulesGC2andGC1.rds")


#### EDIT META FILE

virus <- as.character (gse68310meta$illness)
influenza <- grep ("influenza", virus )
notinfluenza <- setdiff (c (1:880), influenza)

days <- gse68310meta$time
levels (days)
##  "Baseline" "Day0"     "Day2"     "Day4"     "Day6"     "Day21"    "Spring"  
days2 <- days
levels (days2) <- c ("healthy","ill","ill","ill","ill","healthy","healthy")
healthy <- which (days2 == "healthy")
ill <- which (days2 == "ill")

illinfluenza <- intersect (influenza, ill)
illnotinfluenza <- intersect (notinfluenza, ill)
healthyinfluenza <- intersect (influenza, healthy)
healthynotinfluenza <- intersect (notinfluenza, healthy)

illness <- rep ("other", times = 880)
illness [healthyinfluenza] <- "H1"
illness [healthynotinfluenza] <- "H2"
illness [illinfluenza] <- "ill infl"
illness [illnotinfluenza] <- "ill no infl"

gse68310meta2 <- gse68310meta

gse68310meta2$illness <- illness

#### NAME META FILE MORE GENERALLY

metadfnr2 <- gse68310meta2

#### EDIT EXPRESSION TABLE (put sample order as in meta file)

gse68310$ID_REF <- as.character (gse68310$ID_REF)
gse68310$IDENTIFIER <- as.character (gse68310$IDENTIFIER)  
samplesmeta2 <-as.character (metadfnr2$sample)
samplestable <- colnames (gse68310) [3:882]
sampleorder <- numeric (); for (i in samplesmeta2){thisone <- which (samplestable == i); sampleorder <- c (sampleorder, thisone)}
gse68310nr2 <- gse68310 [, c (1,2,sampleorder +2)]


#### NAME EXPRESSION TABLE MORE GENERALLY

tabledfnr2 <- gse68310nr2
nrofcolumns <- dim (tabledfnr2) [2]


####  FOR PLOTTING EXPRESSION VALUES OF GENE SETS GC-2 AND GC-1 (PLOT 1)


##  expression table from log2 scale to linear scale

funexp <- function (x) {
  2 ^ x
}

tablenr2exp <- apply (tabledfnr2  [,3:nrofcolumns],2,funexp)

### subtract too high BG value in expression table

backgroundremove <- rep (150, times = 880)
tablenr2exp <-  tablenr2exp -backgroundremove


tablesize <- dim (tabledfnr2)
samplenr <- tablesize [2] - 2

modulelist <- modules$modulegenesgpl10558idrefs
modulenames <- names (modulelist)
modulenr <- length (modulelist)

df_forplot1 <- as.data.frame (1:samplenr)

counter <- 0


for (k in 1:modulenr) {
  counter <- counter + 1 ; moduleidrefs <- modulelist [[k]]
  
  ## use gene idrefs present on platform
  
  
  presentidrefs <- intersect (moduleidrefs, tabledfnr2$ID_REF)
  these <-
    numeric ();for (xx in presentidrefs) {
      this <- which (tabledfnr2$ID_REF == xx); these <- c (these, this)
    }
  
  
  ##  put expression range for each gene i in module k in a variable named ranges
  
  ranges <- numeric () ; for (i in presentidrefs)  {
    thisone <- which (tabledfnr2$ID_REF == i)
    range <-
      max (tablenr2exp [thisone,]) - min ((tablenr2exp [thisone,]))
    ranges <- c (ranges, range)
  }
  
  
  ## in case a gene module contains only one gene
  
  if (length(these) == 1) {
    these <- c (these,these)
  }
  
  ##  dataframe to collect module info
  
  moduleinfodf <-
    as.data.frame(cbind (
      as.character (tabledfnr2$IDENTIFIER [these]), as.character (tabledfnr2$ID_REF[these]), ranges
    ))
  colnames (moduleinfodf) [1:2] <- c ("gene","ID-REF")
  
  ##  take weights for each gene in module from expression ranges
  
  fractionsofone <- ranges / sum (ranges)
  weights <- 1 / fractionsofone
  moduleinfodf$weights <- weights
  
  
  
  moduleinfodf$genenrs <- these
  
  
  ## use a table for weighted expression values , divide by number of genes in module
  ## number of genes in a module is a constant, used to keep module values more in range of single gene
  ## expression values (dividing by nr of genes in module can be omitted)
  
  weightedtable <-
    tablenr2exp [moduleinfodf$genenrs,] * moduleinfodf$weights / (length (presentidrefs))
  
  
  rownrWT <- dim (weightedtable) [1]
  
  ##  subtract minimal basal expression and platform background for each gene here in case of microarray data
  
  minimalvalues <-
    numeric (); for (p in 1:rownrWT) {
      minimalvalue <-
        min (weightedtable[p,],na.rm = TRUE); minimalvalues <-
          c (minimalvalues,minimalvalue)
    }
  weightedtable <- weightedtable - minimalvalues
  
  
  ##### take module k expression values
  
  ## divide again by number of genes in module (colMeans instead of ColSums), a constant to keep module values more in range of single gene
  ## expression values (can be omitted, then use ColSums)
  
  modulevalues <- colMeans (weightedtable,na.rm = TRUE)
  
  
  
  
  df_forplot1$modulevalues <- as.numeric (as.character (modulevalues))
  colnames (df_forplot1) [counter + 1] <- modulenames [counter]
}


colnames (df_forplot1) [1] <- "sample"

##  add meta info on illness, order of samples in tabledfnr2 and metafile metadfnr2 were made identical


df_forplot1$illness <- metadfnr2$illness
df_forplot1 <- df_forplot1 [,-1]
rownames (df_forplot1 ) <- samplesmeta2


####  FOR PLOTTING EXPRESSION VALUES OF MARKER GENES (PLOT 2)

genes <- c ("LEF1","GYG1","FAM20A","GBP1","OASL")
illuminaidrefs <- c ("ILMN_2213136","ILMN_2230862","ILMN_1812091","ILMN_2148785","ILMN_1674811")

keepnumbers <- numeric (); for (i in illuminaidrefs){thisone <- which (tabledfnr2$ID_REF == i); keepnumbers<- c (keepnumbers, thisone)}
tabledfnr3  <- tabledfnr2 [keepnumbers,]
tablenr3exp <- apply (tabledfnr3  [,3:nrofcolumns],2,funexp)
backgroundremove <- rep (150, times = 880)
tablenr3exp <-  tablenr3exp -backgroundremove
ttablenr3exp <- t (tablenr3exp)
colnames (ttablenr3exp) <- tabledfnr3$IDENTIFIER
ttablenr3expdf <-as.data.frame (ttablenr3exp)
df_forplot2 <- ttablenr3expdf 
df_forplot2$illness <- metadfnr2$illness

#### FOR PLOTTING EXPRESSION VALUES OF GENES IN SET GC-2 (PLOT 3)

genes <- modules$modulegenes$GC_2
illuminaidrefs <- modules$modulegenesgpl10558idrefs$GC_2


keepnumbers <- numeric (); for (i in illuminaidrefs){thisone <- which (tabledfnr2$ID_REF == i); keepnumbers<- c (keepnumbers, thisone)}
tabledfnr3  <- tabledfnr2 [keepnumbers,]
tablenr3exp <- apply (tabledfnr3  [,3:nrofcolumns],2,funexp)
backgroundremove <- rep (150, times = 880)
tablenr3exp <-  tablenr3exp -backgroundremove
ttablenr3exp <- t (tablenr3exp)
colnames (ttablenr3exp) <- tabledfnr3$IDENTIFIER
ttablenr3expdf <-as.data.frame (ttablenr3exp)
df_forplot3 <- ttablenr3expdf 
df_forplot3$illness <- metadfnr2$illness

################# CHOOSE GENES AND GENE SETS FOR PLOTTING ###################

## EITHER marker genes and genesets GC-1 and GC-2 OR single genes in set GC-2

df_forplot <- cbind (df_forplot2 [, c (1,2,3)], df_forplot1 [,c (1,2)], df_forplot2 [, c (4,5,6)])
## df_forplot <- df_forplot3

################# REMOVE ILL NO INFLUENZA AND CORRESPONDING HEALTHY CONTROLS 

healthycontrolsforinfluenza <-  which (df_forplot$illness == "H1")
df_forplot$illness [ healthycontrolsforinfluenza] <- "H"
withinfluenza <-  which (df_forplot$illness == "ill infl")
df_forplot$illness [ withinfluenza] <- "influenza"
remove1 <- which (df_forplot$illness == "H2")
remove2 <- which (df_forplot$illness == "ill no infl")
df_forplot <- df_forplot [-c (remove1, remove2),]


################# TO PLOT ###################################################

plotRstudio_list <- list ()
plotPDF_list <- list ()

nrofmodules <- dim (df_forplot) [2] - 1

modulesforplot <- colnames (df_forplot) [1:nrofmodules]

## plotting each module expression distribution at log2 scale with same Y axis interval of 9 units
## useful for comparing results for multiple datasets with platform GPL10558

for (i in 1:nrofmodules) {
  values  <- df_forplot [,i]
  log2values <- log2 (values)
  log2valuesrange <- max (log2values) - min (log2values)
  halfrange <- log2valuesrange / 2
  midpoint <- min (log2values) +  halfrange
  
  maxpoint <- midpoint + 4.5
  minpoint <- midpoint - 4.5
  titlename <- modulesforplot [i]
  
  ##  if larger Y axis interval used, print gene name with exclamation mark
  
  if (halfrange > 4.5) {
    maxpoint <-
      midpoint + halfrange; minpoint <-
        midpoint - halfrange; titlename <-
          paste (modulesforplot [i], "!" , sep = " ")
  }
  
  
  
  p1 <- "nothing"
  p2 <- "nothing"
  
  ##  use aes_string not aes for plots in list
  
  
  p1 =  ggplot (df_forplot,  aes_string(y = log2values, x = "illness")) +
    
    geom_boxplot(
      width = 0.5,outlier.shape = NA,lwd = 1.5, fatten = 1.5, color = "green"
    ) +
    geom_dotplot(
      binaxis = "y", stackdir = "center", binwidth = 0.1 , alpha = 0.5, color = NA
    ) +
    scale_y_continuous(limits = c (minpoint,maxpoint)) +
    ggtitle(titlename) +
    theme (
      axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size =
                                                                                                  20)
    )
  
  p2  =  ggplot (df_forplot,  aes_string(y = log2values, x = "illness")) +
    
    geom_boxplot(
      width = 0.5,outlier.shape = NA,lwd = 1.5, fatten = 3, color = "green"
    ) +
    geom_dotplot(
      binaxis = "y", stackdir = "center", binwidth = 0.1 , alpha = 0.5, color = NA
    ) +
    scale_y_continuous(limits = c (minpoint,maxpoint)) +
    ggtitle(titlename) +
    theme (
      plot.margin =  unit (c (1,1,1,1),"cm"),plot.title = element_text(face = "bold", size = 80),axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(face = "bold", size =
                                                                                                                                                                                             80), axis.text.y = element_text(face = "bold", size = 40)
    )
  
  
  
  
  
  plotRstudio_list [[i]] <- p1
  plotPDF_list  [[i]] <- p2
  
  
  
}

##  view in RStudio

grid.arrange (grobs = plotRstudio_list)



##  for saving corresponding PDF file
##  3 groups, width is nrofmodules*10


setwd("~/GitHub/GCblood_repo/plotting/plotpdfs")

filenamesindir <- dir ()
nextnr <- length (filenamesindir) + 1
newfilename <- paste ("new", nextnr,"gse68310", ".pdf", sep = "")

pdf(newfilename, width = nrofmodules * 10, height = 12)

grid.arrange (grobs = plotPDF_list, nrow = 1)

dev.off ()

