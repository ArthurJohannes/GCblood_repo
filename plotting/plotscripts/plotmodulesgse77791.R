##### PLOT GENE MODULES EXPRESSION DISTRIBUTIONS FOR MICROARRAY GSE77791 BURN SHOCK

library (gridExtra)
library (ggplot2)


##  read gene expression table, metadata, and gene modules to use

setwd("~/GitHub/GCblood_repo/plotting/plotdata")
gse77791 <- readRDS ("Tablegse77791")
gse77791meta <- readRDS ("gse77791meta.rds")
modules <- readRDS ("modules.rds")

##  EDITING meta file

### remove 168h samples from meta file
 
 samplesnotused <- which (gse77791meta$time == "168h after treatment")
 gse77791meta2 <- gse77791meta [-samplesnotused,] 
 dim (gse77791meta2)
 rownames (gse77791meta2) <- 1:99
 
  print (gse77791meta2)
 
### add new factor illness to meta file
  
  
  illness <- rep ("other", times = 99)
  illness [1:13] <- "healthy"
  illness [14:73] <- "b.shock"
  illness [74:99] <- "b.shock GC"
  gse77791meta2$illness <- illness
  gse77791meta2$illness <- factor (gse77791meta2$illness, levels = c ("healthy","b.shock","b.shock GC"))

  ## EDITING table file
  
  ### subset and order samples in expression table according to edited meta file  
  
  samplesmeta2 <-as.character (gse77791meta2$sample)
  samplestable <- colnames (gse77791) [3:119]
  sampleorder <- numeric (); for (i in samplesmeta2){thisone <- which (samplestable == i); sampleorder <- c (sampleorder, thisone)}
  gse77791nr2 <- gse77791 [, c (1,2,sampleorder +2)]
  

  

nrofcolumns <- dim (gse77791nr2) [2]

summary (colMeans (gse77791nr2 [,3:nrofcolumns],na.rm = TRUE))

##  expression table from log2 scale to linear scale

funexp <- function (x) {
  2 ^ x
}

gse77791nr2exp <- apply (gse77791nr2 [,3:nrofcolumns],2,funexp)
tgse77791nr2exp <- t (gse77791nr2exp)
colnames (tgse77791nr2exp) <- as.character (gse77791nr2$IDENTIFIER)


tablesize <- dim (gse77791nr2)

colnames (gse77791nr2) [1:2] <- c ("ID_REF","IDENTIFIER")
samplenr <- tablesize [2] - 2

##  use pre-selected adequate affymetrix GPL570 idrefs, one idref for each gene in module

modulelist <- modules$modulegenesgpl570idrefs
modulenames <- names (modulelist)
modulenr <- length (modulelist)


samplenr <- tablesize [2] - 2

df_forplot <- as.data.frame (1:samplenr)

counter <- 0


for (k in 1:modulenr) {
  counter <- counter + 1 ; moduleidrefs <- modulelist [[k]]
  
  ## use gene idrefs present on platform
  
  
  presentidrefs <- intersect (moduleidrefs, gse77791nr2$ID_REF)
  these <-
    numeric ();for (xx in presentidrefs) {
      this <- which (gse77791nr2$ID_REF == xx); these <- c (these, this)
    }
  
  
  ##  put expression range for each gene i in module k in a variable named ranges
  
  ranges <- numeric () ; for (i in presentidrefs)  {
    thisone <- which (gse77791nr2$ID_REF == i)
    range <-
      max (gse77791nr2exp [thisone,]) - min ((gse77791nr2exp [thisone,]))
    ranges <- c (ranges, range)
  }
  
  
  ## in case a gene module contains only one gene
  
  if (length(these) == 1) {
    these <- c (these,these)
  }
  
  ##  dataframe to collect module info
  
  moduleinfodf <-
    as.data.frame(cbind (
      as.character (gse77791nr2$IDENTIFIER [these]), as.character (gse77791nr2$ID_REF[these]), ranges
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
    gse77791nr2exp [moduleinfodf$genenrs,] * moduleinfodf$weights / (length (presentidrefs))
  
  
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
  
  
  
  
  df_forplot$modulevalues <- as.numeric (as.character (modulevalues))
  colnames (df_forplot) [counter + 1] <- modulenames [counter]
}


colnames (df_forplot) [1] <- "sample"

##  add meta info on illness, order of samples in table gse77791nr2 and metafile gse77791nr2 were made identical


df_forplot$illness <- gse77791meta2$illness
df_forplot <- df_forplot [,-1]



#################  TO PLOT ############

plotRstudio_list <- list ()
plotPDF_list <- list ()

nrofmodules <- dim (df_forplot) [2] - 1

modulesforplot <- colnames (df_forplot) [1:nrofmodules]

## plotting each module expression distribution at log2 scale with same Y axis interval of 9 units
## useful for comparing results for multiple datasets with platform GPL570

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
newfilename <- paste ("new", nextnr,"gse77791", ".pdf", sep = "")

pdf(newfilename, width = nrofmodules * 10, height = 12)

grid.arrange (grobs = plotPDF_list, nrow = 1)

dev.off ()
