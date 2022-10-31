#### PLOTTING GENE EXPRESSION VALUE DISTRIBUTIONS RNAseq gse154918 #####


library (gridExtra)
library (ggplot2)

setwd("~/GitHub/GCblood_repo/plotting/plotdata")
gse154918 <- readRDS ("Tablegse154918")
gse154918meta <- readRDS ("gse154918meta.rds")

nrofcolumns <- dim (gse154918) [2]

summary (colMeans (gse154918 [,3:nrofcolumns],na.rm = TRUE))

##  expression table from log2 to linear scale

funexp <- function (x) {
  2 ^ x
}

gse154918exp <- apply (gse154918 [,3:nrofcolumns],2,funexp)
tgse154918exp <- t (gse154918exp)
colnames (tgse154918exp) <- as.character (gse154918$IDENTIFIER)

##  select genes for plotting expression value distributions

selectedgenes <-
  c (
    "MARCO","CDKN1C","C1QA","ALDH1A1","ASGR2","PPARG","MERTK","FLT3","ADAMTS2","VSIG4","GBP5","OASL","AZU1","MMP8","HP","MME","GYG1","ARG1","OLAH","DAAM2"
  )

allgenesonplatform <- as.character (gse154918$IDENTIFIER)
keep <-
  numeric  (); for (i in selectedgenes) {
    this <- which (allgenesonplatform == i); keep <- c (keep, this)
  }

tgse154918expkeep <- tgse154918exp [,keep]

##  set groups for comparing expression value distributions
##  sample order in expression table is identical to sample order in meta file


levels (gse154918meta$characteristics_ch1)
illness <- gse154918meta$characteristics_ch1

levels (illness) <- c ("healthy","ill","ill","ill","ill","ill")



df_forplot <- as.data.frame (tgse154918expkeep)
df_forplot$illness <- illness


##########################  TO PLOT  ################################################

plotRstudio_list <- list ()
plotPDF_list <- list ()

nrofgenes <- dim (df_forplot) [2] - 1

genesforplot <- colnames (df_forplot) [1:nrofgenes]

## plotting each gene expression distribution at log2 scale with same Y axis interval of 9 units
## useful for comparing results between genes (expression range) and between multiple RNAseq datasets

for (i in 1:nrofgenes) {
  values  <- df_forplot [,i]
  log2values <- log2 (values)
  log2valuesrange <- max (log2values) - min (log2values)
  halfrange <- log2valuesrange / 2
  midpoint <- min (log2values) +  halfrange
  
  maxpoint <- midpoint + 4.5
  minpoint <- midpoint - 4.5
  titlename <- genesforplot [i]
  
  ##  if larger Y axis interval used, print gene name with exclamation mark
  
  if (halfrange > 4.5) {
    maxpoint <-
      midpoint + halfrange; minpoint <-
        midpoint - halfrange; titlename <-
          paste (genesforplot [i], "!" , sep = " ")
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
newfilename <- paste ("new", nextnr,"gse154918", ".pdf", sep = "")

pdf(newfilename, width = nrofgenes * 6, height = 12)

grid.arrange (grobs = plotPDF_list, nrow = 2)

dev.off ()
