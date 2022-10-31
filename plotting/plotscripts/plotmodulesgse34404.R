library (gridExtra)
library (ggplot2)

##  read gene expression table, metadata, and gene modules to use


setwd("~/GitHub/GCblood_repo/plotting/plotdata")

gse34404 <- readRDS ("Tablegse34404")
gse34404meta <- readRDS ("gse34404meta.rds")
modules <- readRDS ("modules.rds")



nrofcolumns <- dim (gse34404) [2]

summary (colMeans (gse34404 [,3:nrofcolumns],na.rm = TRUE))

##  expression table gse34404 is already linear scale

funidentity <- function (x) {
  x
}

gse34404exp <- apply (gse34404 [,3:nrofcolumns],2,funidentity)
tgse34404exp <- t (gse34404exp)
colnames (tgse34404exp) <- as.character (gse34404$IDENTIFIER)


tablesize <- dim (gse34404)

colnames (gse34404) [1:2] <- c ("ID_REF","IDENTIFIER")
samplenr <- tablesize [2] - 2

##  use pre-selected adequate illumina GPL10558 idrefs, one idref for each gene in module

modulelist <- modules$modulegenesgpl10558idrefs
modulenames <- names (modulelist)
modulenr <- length (modulelist)


samplenr <- tablesize [2] - 2

df_forplot <- as.data.frame (1:samplenr)

counter <- 0


for (k in 1:modulenr) {
  counter <- counter + 1 ; moduleidrefs <- modulelist [[k]]
  
  ## use gene idrefs present on platform
  
  
  presentidrefs <- intersect (moduleidrefs, gse34404$ID_REF)
  these <-
    numeric ();for (xx in presentidrefs) {
      this <- which (gse34404$ID_REF == xx); these <- c (these, this)
    }
  
  
  ##  put expression range for each gene i in module k in a variable named ranges
  
  ranges <- numeric () ; for (i in presentidrefs)  {
    thisone <- which (gse34404$ID_REF == i)
    range <-
      max (gse34404exp [thisone,]) - min ((gse34404exp [thisone,]))
    ranges <- c (ranges, range)
  }
  
  
  ## in case a gene module k contains only one gene
  
  if (length(these) == 1) {
    these <- c (these,these)
  }
  
  ##  dataframe to collect module info
  
  moduleinfodf <-
    as.data.frame(cbind (
      as.character (gse34404$IDENTIFIER [these]), as.character (gse34404$ID_REF[these]), ranges
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
    gse34404exp [moduleinfodf$genenrs,] * moduleinfodf$weights / (length (presentidrefs))
  
  
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

##  add meta info on illness, order of samples in table and metafile gse34404 are identical

levels (gse34404meta$illness)
illness <- gse34404meta$illness

levels (illness) <- c ("control","sympt.")
df_forplot$illness <- illness
df_forplot <- df_forplot [,-1]




#################  TO PLOT ############

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
  
  ##  if larger Y axis interval used, print module name with exclamation mark
  
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

setwd("~/GitHub/GCblood_repo/plotting/plotpdfs")

filenamesindir <- dir ()
nextnr <- length (filenamesindir) + 1
newfilename <- paste ("new", nextnr,"gse34404", ".pdf", sep = "")

pdf(newfilename, width = nrofmodules * 6, height = 12)

grid.arrange (grobs = plotPDF_list, nrow = 1)

dev.off ()
