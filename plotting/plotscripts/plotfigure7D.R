####### PLOT GENE (SETS) EXPRESSION VALUE DISTRIBUTION FOR MICROARRAY GSE42026 SEVERE RESP. VIR. INFECTION ##

library (gridExtra)
library (ggplot2)


##  read gene expression table, metadata, and gene modules to use

setwd("~/GitHub/GCblood_repo/plotting/plotdata")

gse42026 <- readRDS ("Tablegse42026")
gse42026meta <- readRDS ("gse42026meta.rds" )
modules <- readRDS ("modulesGC2andGC1.rds")

#### EDIT META FILE

 levels (gse42026meta$illness)
                              
 levels (gse42026meta$illness) <- c ("bact","infl","H","RSV")
 gse42026meta$illness <- as.character(gse42026meta$illness)
 gse42026meta$illness <- factor(gse42026meta$illness, levels = c ("H","bact","infl","RSV"))
 

gse42026meta2 <- gse42026meta

#### NAME META FILE MORE GENERALLY

metadfnr2 <- gse42026meta2

#### EDIT EXPRESSION TABLE (put sample order as in meta file)

nrofcolumns <- dim (gse42026) [2]
gse42026$ID_REF <- as.character (gse42026$ID_REF)
gse42026$IDENTIFIER <- as.character (gse42026$IDENTIFIER)  

samplesmeta2 <-as.character (metadfnr2$sample)
samplestable <- colnames (gse42026) [3:94]
sampleorder <- numeric (); for (i in samplesmeta2){thisone <- which (samplestable == i); sampleorder <- c (sampleorder, thisone)}
gse42026nr2 <- gse42026 [, c (1,2,sampleorder +2)]

#### NAME EXPRESSION TABLE MORE GENERALLY

tabledfnr2 <- gse42026nr2 

summary (colMeans (tabledfnr2 [,3:nrofcolumns],na.rm = TRUE))
### expression values already linear, add value to have positive values only 

funplus15 <- function (x) {
  x +15
}

tablenr2exp <- apply (tabledfnr2  [,3:nrofcolumns],2,funplus15)

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
tablenr3exp <- tablenr2exp [keepnumbers,]
ttablenr3exp <- t (tablenr3exp)
colnames (ttablenr3exp) <-  tabledfnr2$IDENTIFIER [keepnumbers]
ttablenr3expdf <-as.data.frame (ttablenr3exp)
df_forplot2 <- ttablenr3expdf 
df_forplot2$illness <- metadfnr2$illness

#### FOR PLOTTING EXPRESSION VALUES OF GENES IN SET GC-2 (PLOT 3)

genes <- modules$modulegenes$GC_2
illuminaidrefs <- modules$modulegenesgpl10558idrefs$GC_2

keepnumbers <- numeric (); for (i in illuminaidrefs){thisone <- which (tabledfnr2$ID_REF == i); keepnumbers<- c (keepnumbers, thisone)}
tablenr3exp <- tablenr2exp [keepnumbers,]
ttablenr3exp <- t (tablenr3exp)
colnames (ttablenr3exp) <- tabledfnr2$IDENTIFIER [keepnumbers]
ttablenr3expdf <-as.data.frame (ttablenr3exp)
df_forplot3 <- ttablenr3expdf 
df_forplot3$illness <- metadfnr2$illness

################# CHOOSE GENES AND GENE SETS FOR PLOTTING ###################

## EITHER marker genes and genesets GC-1 and GC-2 OR single genes in set GC-2

df_forplot <- cbind (df_forplot2 [, c (1,2,3)], df_forplot1 [,c (1,2)], df_forplot2 [, c (4,5,6)])
## df_forplot <- df_forplot3


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
      axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size =
                                                                                                  40), axis.text.y = element_text(size = 40),plot.title = element_text(size =
                                                                                                                                                                         60)
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
newfilename <- paste ("new", nextnr,"gse42026", ".pdf", sep = "")

pdf(newfilename, width = nrofmodules * 10, height = 12)

grid.arrange (grobs = plotPDF_list, nrow = 1)

dev.off ()


