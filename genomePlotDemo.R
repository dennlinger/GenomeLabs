# genomePlotDemo.R
#
# Purpose:  Demo a genome plot of genes in a circle plot with functional
#           connections.
#
#
# Version:  1.0
# Date:     2018 03 17
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Dependencies:
#           readr package
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Version history:
#   1.0  Final version for Biohacks 2018
#   0.5  Added an Introduction section and more comments
#   0.4  Update for updated datafiles
#   0.3  Bugfix in coordinate scaling to SVG coordinates for rectangles
#   0.2  Improve abstractions and modularization, move functions to
#        separate file.
#   0.1  First draft
#
# ToDo:
#    - ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                   Line
#TOC> -------------------------------------------------------
#TOC>   1        INTRODUCTION                              55
#TOC>   2        PARAMETERS                                93
#TOC>   3        PACKAGES AND FUNCTIONS                   118
#TOC>   4        PROCESS                                  127
#TOC>   4.1        READ SOURCE DATA                       133
#TOC>   4.2        INITIALIZE DATA STRUCTURES             142
#TOC>   4.3        ANNOTATE                               181
#TOC>   4.3.1          Annotate relationship types:       193
#TOC>   4.3.2          Annotate relationship weights:     210
#TOC>   4.4        LAYOUT                                 230
#TOC>   4.5        PLOT                                   377
#TOC>   4.5.1          Compute scale and translation      380
#TOC>   4.5.2          Write SVG header                   426
#TOC>   4.5.3          Render all elements                431
#TOC>   4.5.4          Write SVG footer                   442
#TOC>   5        FINISH                                   446
#TOC> 
#TOC> ==========================================================================


# =    1  INTRODUCTION  ========================================================

# You should familiarize yourself with genomePlotBasic.R before you work though
# this script - this script provides an extension to the basic workflow. Genes
# are plotted on a circular chromosome, and edges are added between genes that
# share a GO annotation from the biological process GOslim ontology. Further
# extensions of the code are provided in genomePlotIntermediate.R (There is no
#  "advanced" version - that would be your code).

# The code proceeds through five steps:

#  1 -  Read the source data:
#         Just like in genomePlotBasic.R a single data file containing gene
#         annotations is read into a data frame.
#
#  2 -  Initialize data structures:
#         In this step, information objects are defined. In this demo we define
#         gene information in a data frame, just like we did in
#         genomePlotBasic.R.  Then we define edges between
#         genes that are annoatetd to the same GO term.
#
#  3 -  Annotate:
#         genomePlotBasic.R used only data that was read from the source file.
#         In this demo, we compute annotations for relationships. Each
#         relationship gets a category, and a weight.
#
#  4 -  Layout:
#         genomePlotBasic.R mapped genes to a linear chromosome. This script
#         plots each gene as a coloured rectangle on a circle, and draws lines
#         for each relationship that we have defined above.
#
#  5 -  Plot:
#         The plotting step is very similar to genomePlotBasic.R, except for
#         a bit of code to collect the boundaries of the shapes we will draw
#         in order to scale our plot into the available window.



# =    2  PARAMETERS  ==========================================================

# Code should not contain "magic numbers". Constants that we need are
# defined and commented here.

CHR20LENGTH   <- 64444167  # basepairs of the chromosome we are working with

DATAFILE <- "Chr20GeneData.tsv"  # Chromosome data input. See README-DATA and
                                 # prepareGenomeData.R for a description of the
                                 # contents, and the code that produces it from
                                 # database sources.

NACOLOUR <- "#AAAAAA"            # Neutral grey for NA attributes

SVGFILE <- "interactions.svg"  # Filename for the output we produce

# UTPoster prints from 24" x 36" all the way to 60" x 300".
# Let's assume letter size for this demo, and subtract a 1" margin on both
# sides.
PAGEWIDTH  <- ( 8.5 - 2) * 2.54    # in cm
PAGEHEIGHT <- (11.0 - 2) * 2.54    # in cm
RESOLUTION <- 150                  # pixels per 2.54 cm



# =    3  PACKAGES AND FUNCTIONS  ==============================================
#
# All required packages and functions are loaded from the source file below.
# You can inspect/copy/modify the source code there.
setwd("/home/dennis/BioHacks/2018-Challenge/")
source("genomePlotFunctions.R")



# =    4  PROCESS  =============================================================

# This demo code will plot genes as rectangles on a circle, color the boxes, and
# connect genes that share the same function category with a line.


# ==   4.1  READ SOURCE DATA  ==================================================

# read_tsv() is from the readr package. It is similar to base R's read.delim()
# function, but more modern. Note that it returns a "tibble" which is similar
# but not identical to a data frme.

myData <- read_tsv(DATAFILE)


# ==   4.2  INITIALIZE DATA STRUCTURES  ========================================

# There are many possibilities to store the data for the objects we will analyze
# and draw. Here we take a very simple approach and store gene-level data in one
# data frame, relationship data in another data frame. Gene level data is
# populated directly from the source data, relationship data (edges) are
# computed from the gene data.

# Entities and attributes: data for each gene
myGenes <- data.frame(sym = myData$sym,            # Gene symbols
                      start = myData$start,        # start
                      end = myData$end,            # end
                      strand = myData$strand,      # strand
                      GOid = myData$GO_P,          # GO annotation for "Process"
                      stringsAsFactors = FALSE)


# Relationship annotations: define an edge from <symbol> to <symbol>
myEdges <- data.frame(from = character(),
                      to = character(),
                      stringsAsFactors = FALSE)

for (i in 1:nrow(myGenes)){
  # for each gene, define an edge to all other genes with the same GO ID.
  thisSym <- myGenes$sym[i]
  thisGOid <- myGenes$GOid[i]

  if (! is.na(thisGOid)) {
    sel <- which(myGenes$GOid == thisGOid) # all genes with this GO id
    sel <- sel[sel != i] # remove the index of the original (no self-edge)

    myEdges <- rbind(myEdges, data.frame(from = rep(thisSym, length(sel)),
                                         to = myGenes$sym[sel],
                                         stringsAsFactors = FALSE))
  }
}



# ==   4.3  ANNOTATE  ==========================================================

# We will derive the following annotations from the data we have loaded:
#
# A: each relationship will get a "type". For this demo the type is simply
#    the same as the GO annotation. It could be anything though.
# B: each relationship will get a "weight". For this demo the weight is simply
#    1 minus the frequency of the GO annotation, that is: more specific
#    annotations get a lower weight, more generic annotations get a higher
#    weight.


# ===   4.3.1  Annotate relationship types:  

myEdges$type <- character(nrow(myEdges)) # Add a "type" column - this could
                                         # be any kind of categorical data,
                                         # in this demo we simply use the GO ID.

for (i in 1:nrow(myEdges)){  # for each edge, add the GO id of the "from" gene.
  # Note: we don't need to test for NA here, because we only used annotated
  # genes to build the list.

  sel <- which(myGenes$sym == myEdges$from[i])[1]  # take first match
  myEdges$type[i] <- myGenes$GOid[sel]
}




# ===   4.3.2  Annotate relationship weights:

# We can encode information in the edge thickness. In this demo, we draw
# more specific (less frequently annotated) GO terms with a thicker line.


wGOA <- linMap(table(myGenes$GOid),
               low = 0.99,
               high = 0.01) # map frequencies to weights

myEdges$weight <- numeric(nrow(myEdges)) # add "weight" column

for (i in 1:nrow(myEdges)){  # for each edge, add the weight of the "from" gene.
  myEdges$weight[i] <- wGOA[myEdges$type[i]]
}


# Done with annotations


# ==   4.4  LAYOUT  ============================================================

# The layout phase is where we turn data into visuals.

# Since we need to accommodate quite different types of objects, we will collect
# them in a list. Each element is itself a list that describes the
# object with all the detail we need so we can draw it out later.

myShapes <- list()

# We will call the things that we are going to draw "shapes".
# To draw our shapes, we need to define:
#  - what they are
#  - where they are going to be drawn
#  - how large
#  - with what stroke-width
#  - with what colour
#  - with what fill
#  - ...  Many other attributes can be added - perspective, labels,
#         curvature, line type, shadow, gradient etc. etc.


# Let's start with some generic elements to structure our plot. We will draw the
# chromosome backbone as a circle and we will add the chromosome name as a text
# element.

# At first, we are only concerned with relative positions and we will layout
# shapes into an arbitrary canvas. Later we will map this into page coordinates.
# The chromosome circle we will define will be centred on (1, 1), and it we will
# give it a radius of 1.0

# Chromosome backbone:

CHR20ORI <- c(1.0, 1.0)
CHR20RAD <- 1.0

myShapes[[1]] <- list(type = "circle",
                      centre = CHR20ORI,
                      radius = CHR20RAD,
                      fill = "#FFFFFF",    # fill colour
                      stroke = "#4499AA",  # colour of outline
                      sw = 7.0)            # stroke-width


# Next we add some descriptive text:

myShapes[[2]] <- list(type = "text",
                      text = "CHR 20",
                      centre = CHR20ORI,
                      size = 48,           # points
                      font = "Times",
                      fill = "#33AAFF")

# Next we add the genes. We will draw them as rectangles, with a height of a
# fraction of the circle radius. We will place them on the circle at their
# fractional position on the chromosome, and we will rotate them so they point
# radially on the origin. We will also give them colour:

# Colour:
# We are providing a function category2colour() to make life simple. It
# takes a vector of items, and returns a named vector of corresponding
# color values.
#
# For this demo we will color the genes by their function. Thus we define
# a basic, divergent spectrum, and map these colors to GO ids.

mySpect <- c("#f2003c",  # red
             "#F0A200",  # orange
             "#f0ea00",  # yellow
             "#62C923",  # green
             "#0A9A9B",  # blue
             "#1958C3",  # indigo
             "#8000D3",  # violet
             "#D0007F")  # red

myGOcolours <- category2colour(sort(unique(myGenes$GOid)),
                               col = mySpect)


# Add each gene to the list:

for (i in 1:nrow(myGenes)) {

  # The centre of the rectangle is placed on the circle
  # coord2circle() returns x, y, and rotation angle
  circDat <- coord2circle(mean(c(myGenes$start[i], myGenes$end[i])),
                          CHR20LENGTH,
                          CHR20ORI,
                          CHR20RAD)
  width <- abs(myGenes$start[i] - myGenes$end[i]) / CHR20LENGTH  # relative ...
  width <- width * 2 * pi * CHR20RAD                             # on the circle

  height <- 0.05 * CHR20RAD   # 5% of circle radius

  # Layout for genes: each gene will be drawn as a colored box, placed on the
  # circle and rotated appropriately. We store the (x, y) of the centre, as well
  # as width and height of the rectangle, and the angle through which it should
  # be rotated.

  # We use the colour for GO anotations we defined above.
  myFill <- myGOcolours[myGenes$GOid[i]]

  if (is.na(myFill)) {
    myFill <- NACOLOUR
  }

  # At this scale, a typical gene is about a hair's width. We draw the outline
  # of the rectangle, with a thin line, to give it a minimum width for
  # visibility.

  # Define a rectangle shape with these parameters:
  myShapes[[length(myShapes) + 1]] <- list(type = "rect",
                                           centre = circDat[1:2],
                                           w = width,
                                           h = height,
                                           ang = circDat[3],
                                           fill = myFill,
                                           stroke = myFill,
                                           sw = 0.5)  # points
}



# Next, add each relationship to the list:

for (i in 1:nrow(myEdges)) {
  iFrom <- which(myGenes$sym == myEdges$from[i])
  iTo <-   which(myGenes$sym == myEdges$to[i])
  xyFrom <- coord2circle(mean(c(myGenes$start[iFrom], myGenes$end[iFrom])),
                         CHR20LENGTH,
                         CHR20ORI,
                         CHR20RAD)[1:2]
  xyTo   <- coord2circle(mean(c(myGenes$start[iTo], myGenes$end[iTo])),
                         CHR20LENGTH,
                         CHR20ORI,
                         CHR20RAD)[1:2]
  myShapes[[length(myShapes) + 1]] <- list(type = "line",
                                           p1 = xyFrom,
                                           p2 = xyTo,
                                           stroke = myGOcolours[myEdges$type[i]],
                                           sw = myEdges$weight[i])
}


# Done. All shapes are defined.


# ==   4.5  PLOT  ==============================================================
# cf. https://www.w3.org/TR/SVG

# ===   4.5.1  Compute scale and translation 

# Caution: the SVG coordinate system has its origin (0, 0) in the TOP LEFT
# corner, positive X goes right, and positive Y goes down. Here we define the
# necessary scaling and translation.

# First: we fetch the centres of genes from the list shapes to compute
# the range of x and y values we will plot.

xs <- numeric()
ys <- numeric()

for (i in 1:length(myShapes)) {
  if (myShapes[[i]]$type == "rect") {
    xs <- c(xs, myShapes[[i]]$centre[1])
    ys <- c(ys, myShapes[[i]]$centre[2])
  }
}


# Next, we compute the range of x and y values:

dX <- range(xs)[2] - range(xs)[1]
dY <- range(ys)[2] - range(ys)[1]


# Given the range that needs to fit on the page, we can compute the scale:

sXY <- min((RESOLUTION * (PAGEWIDTH / 2.54)) / dX,
         (RESOLUTION * (PAGEHEIGHT / 2.54)) / dY)

sXY <- sXY * 0.95 # tweak it a bit smaller to allow for stroke widths


# We compute the dimensions of the page in pixels ...

Xpx <- RESOLUTION * (PAGEWIDTH  / 2.54)
Ypx <- RESOLUTION * (PAGEHEIGHT / 2.54)

# And we compute a translation: for this demo, we move our CHR20ORI
# to the centre of the page:

tXY <- c(Xpx / 2, Ypx / 2)  - (sXY * CHR20ORI)  # translate



# ===   4.5.2  Write SVG header              
mySVG <- SVGheader()
mySVG <- c(mySVG, SVGdefinePage(Xpx, Ypx))


# ===   4.5.3  Render all elements           
#
for (i in 1:length(myShapes)) {

  mySVG <- c(mySVG, SVGrenderElement(myShapes[[i]],
                                     sc = sXY,
                                     tr = tXY,
                                     Y = Ypx))
}


# ===   4.5.4  Write SVG footer              
mySVG <- c(mySVG, SVGfooter())


# =    5  FINISH  ==============================================================

# Write the SVG to file
writeLines(mySVG, con = SVGFILE)

# Open the SVG in the default browser to visualize
system(sprintf("open -a \"Google Chrome\" %s", SVGFILE))   # For MacOS
# Windows ???
# Linux ???



# ====  TESTS  =================================================================
# ...





# [END]
