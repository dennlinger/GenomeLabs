# genomePlotIntermediate.R
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
#   0.2  Implemented Force Directed Layout, and cubic Bezier curve edges
#   0.1  Derived from genomePlotDemo.R   V 0.3
#
# ToDo:
#    - ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                             Line
#TOC> -----------------------------------------------------------------
#TOC>   1        INTRODUCTION                                        58
#TOC>   2        PARAMETERS                                          93
#TOC>   3        PACKAGES AND FUNCTIONS                             124
#TOC>   4        READ SOURCE DATA                                   155
#TOC>   5        INITIALIZE DATA STRUCTURES                         166
#TOC>   6        ANNOTATE                                           201
#TOC>   6.1        Annotate Gene-Protein edge types:                213
#TOC>   6.2        Annotate Gene-Protein edge weights:              225
#TOC>   7        LAYOUT                                             250
#TOC>   7.1        Chromosome backbone                              286
#TOC>   7.2        Descriptive text                                 302
#TOC>   7.3        Layout Genes                                     313
#TOC>   7.4        Layout Proteins                                  387
#TOC>   7.4.1          Force directed layout: setup                 440
#TOC>   7.4.2          Force directed layout: iterations            534
#TOC>   7.4.3          Force directed layout: coordinate update     620
#TOC>   7.5        Layout functional interaction edges              663
#TOC>   8        PLOT                                               706
#TOC>   8.1        Compute scale and translation                    709
#TOC>   8.2        Write SVG header                                 755
#TOC>   8.3        Render all elements                              760
#TOC>   8.4        Write SVG footer                                 771
#TOC>   9        FINISH                                             775
#TOC>
#TOC> ==========================================================================


# =    1  INTRODUCTION  ========================================================

# You should familiarize yourself with genomePlotDemo.R before you work though
# this script - this script provides an extension to the workflow implemented there. Genes
# are plotted on a circular chromosome. Then selected proteins are plotted onto a concentric circle. Finally the layout of proteins is optimized with a force directed layout algorithm, and edges between the proteins are drawn as a cubic bezier curve.

# The code proceeds through five steps:

#  1 -  Read the source data:
#         Just like in genomePlotDemo.R a single data file containing gene
#         annotations is read into a data frame. In addition, protein network
#         data is read.
#
#  2 -  Initialize data structures:
#         In this step, information objects are defined. In this demo we define
#         gene information in a data frame, just like we did in
#         genomePlotDemo.R.  Then we define edges between
#         genes and those proteins for which we have network information.
#         Finally we add edges that correspond to functional interactions.
#
#  3 -  Annotate:
#         Similar to genomePlotDemo.R we compute categories, and weight
#         annotations for edges between genes and proteins.
#
#  4 -  Layout:
#         genomePlotDemo.R mapped genes to a circle. This script maps proteins
#         to a concentric circle, and optimizes their position with
#         force-directed layout. Rather than drawing simple lines, we draw
#         cubic bezier curves to connect proteins that have a functional
#         interaction.
#
#  5 -  Plot:
#         The plotting step is very similar to genomePlotDemo.R.


# =    2  PARAMETERS  ==========================================================

# Code should not contain "magic numbers". Constants that we need are
# defined and commented here.

CHR20LENGTH   <- 64444167  # basepairs of the chromosome we are working with

DATAFILE <- "Chr20GeneData.tsv"  # Chromosome data input. See README-DATA and
                                 # prepareGenomeData.R for a description of the
                                 # contents, and the code that produces it from
                                 # database sources.

FIFILE <- "Chr20FuncIntx.tsv"    # Functional interaction data input. See
                                 # README-DATA for details.

GOFILE <- "Chr20GOslimData.tsv"  # GO term data input. See
                                 # README-DATA for details.

NACOLOUR <- "#AAAAAA"            # Neutral grey for NA attributes

SVGFILE <- "test.svg"  # Filename for the output we produce

# UTPoster prints from 24" x 36" all the way to 60" x 300".
# Let's assume letter size for this demo, and subtract a 1" margin on both
# sides.
PAGEWIDTH  <- ( 8.5 - 2) * 2.54    # in cm
PAGEHEIGHT <- (11.0 - 2) * 2.54    # in cm
RESOLUTION <- 150                  # pixels per 2.54 cm



# =    3  PACKAGES AND FUNCTIONS  ==============================================
#
# Most required packages and functions are loaded from the source file below.
# You can inspect/copy/modify the source code there.

source("genomePlotFunctions.R")

# Force-directed layout can run for a while. A progress bar tells us
# how we are doing ...
pBar <- function(i, l, nCh = 50) {
  # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)[1]  # use only first, in case there are ties
      p1 <- paste(rep("#", p), collapse = "")
      p2 <- paste(rep("-", nCh - p), collapse = "")
      cat(sprintf("\r|%s%s|", p1, p2))
      flush.console()
    }
  }
  else { # done
    cat("\n")
  }
}



# =    4  READ SOURCE DATA  ====================================================

# read_tsv() is from the readr package. It is similar to base R's read.delim()
# function, but more modern. Note that it returns a "tibble" which is similar
# but not identical to a data frme.

myData <- read_tsv(DATAFILE)   # gene data
myFI <- read_tsv(FIFILE)       # functional interactions
myGO <- read_tsv(GOFILE)       # GO terms


# =    5  INITIALIZE DATA STRUCTURES  ==========================================

# There are many possibilities to store the data for the objects we will analyze
# and draw. Here we take a very simple approach and store entities in data
# frames of attributes
#   - genes
#   - proteins
#   - gene-protein edges
#   - protein-protein edges

# Entities and attributes: data for each gene
myGenes <- data.frame(sym = myData$sym,            # Gene symbols
                      start = myData$start,        # start
                      end = myData$end,            # end
                      strand = myData$strand,      # strand
                      GOid = myData$GO_P,          # GO annotation for "Process"
                      stringsAsFactors = FALSE)

# Proteins: taken from unique gene names for functional interactions
myProteins <- data.frame(sym = unique(c(myFI$a, myFI$b)),   # Gene symbols
                         stringsAsFactors = FALSE)


# Gene protein edges: from <symbol> to <symbol>
myGPedges <- data.frame(g = myProteins$sym,
                        p = myProteins$sym,   # They are the same symbol
                        stringsAsFactors = FALSE)


# Functional interaction edges: between <symbol> and <symbol>
myFIedges <- data.frame(a = myFI$a,
                        b = myFI$b,
                        stringsAsFactors = FALSE)


# =    6  ANNOTATE  ============================================================

# We will derive the following annotations from the data we have loaded:
#
# A: each gene-protein edge will get a "type". For this demo the type is simply
#    the same as the GO annotation. It could be anything though.
# B: each gene-protein edge will get a "weight". For this demo the weight is
#    simply 1 minus the frequency of the GO annotation, that is: more specific
#    annotations get a lower weight, more generic annotations get a higher
#    weight.


# ==   6.1  Annotate Gene-Protein edge types:  =================================

myGPedges$type <- character(nrow(myGPedges)) # Add a "type" column - here
                                             # we simply use the GO ID.

for (i in 1:nrow(myGPedges)){  # for each edge, add the GO id of the gene.

  sel <- which(myGenes$sym == myGPedges$g[i])
  myGPedges$type[i] <- myGenes$GOid[sel]
}


# ==   6.2  Annotate Gene-Protein edge weights:  ===============================

# We can encode information in the edge thickness. Here we draw
# weights according to the Shannon information of GO terms, calculated
# from the counts of annotations in the GO data. We only consider GO terms
# from the biological process ("P") ontology.

myGO <- as.data.frame(myGO[myGO$namespace == "P", ])
rownames(myGO) <- myGO$ID

# GO:0008150 is the root term ID of the "biological process" ontology
nAll <- myGO["GO:0008150", "counts"]

# Compute Shannon information as -log2(counts/nAll) for each edge
myGPedges$weight <- -log2(myGO[myGPedges$type, "counts"]/nAll)

# Linear scale information to stroke width (0, 1)
myGPedges$weight <- linMap(myGPedges$weight, low = 0.01, high = 0.99)




# Done with annotations


# =    7  LAYOUT  ==============================================================

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

# Since we will update some shape coordinates with force-directed layout later,
# we give each shape a "category": "FDLignore", for layout elements; "FDLfixed",
# for elements that contribute to force calculations bot won't themselves move;
# "FDLmovable", for elements whose positions we'll update.

# ==   7.1  Chromosome backbone  ===============================================


CHR20ORI <- c(1.0, 1.0)
CHR20RAD <- 1.0

myShapes[[1]] <- list(type = "circle",
                      class = "layout",
                      category = "FDLignore",
                      centre = CHR20ORI,
                      radius = CHR20RAD,
                      fill = "#FFFFFF",    # fill colour
                      stroke = "#4499AA",  # colour of outline
                      sw = 7.0)            # stroke-width


# ==   7.2  Descriptive text  ==================================================

myShapes[[2]] <- list(type = "text",
                      class = "layout",
                      category = "FDLignore",
                      text = "CHR 20",
                      centre = CHR20ORI,
                      size = 48,           # in points
                      font = "Times",
                      fill = "#33AAFF")

# ==   7.3  Layout Genes  ======================================================

# Next we add the genes. We will draw them as rectangles, with a height of a
# fraction of the circle radius. We will place them on the circle at their
# fractional position on the chromosome, and we will rotate them so they point
# radially on the origin. They will contribute to FDL calculations, but won't
# move - i.e. their category is FDLfixed. We will also give them colour:

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


# Add each gene to the shape list:

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
                                           class = "gene",
                                           category = "FDLfixed",
                                           name = myGenes$sym[i],
                                           centre = circDat[1:2],
                                           w = width,
                                           h = height,
                                           ang = circDat[3],
                                           fill = myFill,
                                           stroke = myFill,
                                           sw = 0.5)  # points
}

# ==   7.4  Layout Proteins  ===================================================

# Next, we add proteins on a concentric circle with a radius of
# 2/3 the chromosome radius.

PROTSCALE <- 2/3

# The circle itself first ...
myShapes[[length(myShapes) + 1]] <- list(type = "circle",
                                         class = "layout",
                                         category = "FDLignore",
                                         centre = CHR20ORI,
                                         radius = CHR20RAD * PROTSCALE,
                                         fill = "none",    # fill colour
                                         stroke = "#9944AA",  # colour line
                                         sw = 3.0)            # stroke-width


# then the proteins - we will draw them as small circles:

pRad  <- 0.02 * CHR20RAD

for (i in 1:nrow(myProteins)) {

  # We initially place them on the circle, at the same radial angle as the
  # gene. We'll update these coordinates later.

  iGene = which(myGenes$sym == myProteins$sym[i])
  circDat <- coord2circle(mean(c(myGenes$start[iGene], myGenes$end[iGene])),
                          CHR20LENGTH,
                          CHR20ORI,
                          CHR20RAD * PROTSCALE)

  # We use the same colour as the genes.
  myFill <- myGOcolours[myGenes$GOid[iGene]]

  if (is.na(myFill)) {
    myFill <- NACOLOUR
  }

  # Define a circle shape with these parameters:
  myShapes[[length(myShapes) + 1]] <- list(type = "circle",
                                           class = "protein",
                                           category = "FDLmovable",
                                           name = myProteins$sym[i],
                                           centre = circDat[1:2],
                                           radius = pRad,
                                           fill = myFill,
                                           stroke = "#000000",
                                           sw = 1.0)  # points
}


# ===   7.4.1  Force directed layout: setup

# In this demo code, we will use a force directed layout to arrange the
# proteins on the circle in the following way: we tether the position
# of a protein to its gene. Then we pull proteins with functional interactions
# near each other.

# Forces in such a scheme are typically linear for an attractive term (F = k*d)
# and quadratic for a repulsive term (F = -k/d^2). The attractive term pulls
# connected nodes towards each other, while the repulsive term keeps them from
# overlapping. We calculate the vector between each node, scale it by the
# distance-dependent force, and add all vectors together. Then we move each node
# a little bit along the resultant vector, and iterate.

# There is a slight complication here, since we constrain movement on a circle.
# We thus move our nodes in small increments and project them back on the circle
# after each step. This specific part of the code is easily eliminated for the
# general case. However, nodes can't pass each other on the circle in this
# scheme, and thus we will run the layout computation in two schedules: at first
# switching repulsive forces  off, then gradually switching the repulsive
# forces back on and equilibrating some more. How to schedule the layout well is
# a matter of trial and error and what we are doing in this code is unlikely to
# be the best solution. Our example is just simple to code and - we hope - easy
# to understand.

# Disclaimer: normally we would put much of the code we present below
# into functions. In this demo code, we keep it exposed to make the process
# more explicit.

# Extract the coordinate information from the shape list:
fdl <- data.frame(name = character(),           # the gene name,
                  x    = numeric(),             # its current x ...
                  y    = numeric(),             # ... and y coordinate,
                  isMovable = logical(),        # whether it should move.
                  stringsAsFactors = FALSE)

for (i in 1:length(myShapes)) {

  thisCat <- myShapes[[i]]$category
  if (thisCat %in% c("FDLfixed", "FDLmovable") &&   # use relevant shapes that
      myShapes[[i]]$name %in% myGPedges$g) {        # are also in myGPedges
    idx <- nrow(fdl) + 1
    fdl[idx, "name"] <- myShapes[[i]]$name
    fdl[idx, c("x", "y")] <- myShapes[[i]]$centre
    # proteins are movable, genes are not
    fdl[idx, "isMovable"] <- ifelse(thisCat == "FDLmovable", TRUE, FALSE)
  }
}

# Define force-functions for attractive and repulsive terms. Not all pairs
# of nodes attract each other, but all of them are potentially repulsive.
# In large layout problems, we would limit calculation of repulsive forces
# by coarse-graining (no need to calculate long distance forces for far
# away points on every iteration), and distance cutoff.

fA <- function(p, v, kA) {
  # calculate a force vector on p, exerted from points in the vector v, with
  # a linear attractive term (F = k*d)
  #
  if (is.vector(v)) {
    v <- matrix(v, nrow = 1)
  }
  return(kA * colSums(v-p))   # sum and scale
}

fR <- function(p, v, kR) {
  # calculate a force vector on p, exerted from points in the vector v, with
  # a quadratic repulsive term (F = -k/d^2)

  if (is.vector(v)) {
    v <- matrix(v, nrow = 1)
  }

  v[,1] <- v[,1] - p[1]               # shift to origin
  v[,2] <- v[,2] - p[2]
  d <- sqrt( (v[,1])^2 + (v[,2])^2 )  # distances
  F <- -kR / d^2                      # forces
  v <- v / d                          # normalize
  v <- F * v                          # scale
  return(colSums(v))                  # sum
}


# define a function to project a point on a circle.
toCircle <- function(p, ori, rad) {
  # project point p onto the circle with centre ori and radius rad.
  p <- p - ori
  d <- sqrt((p[1])^2 + (p[2])^2)
  p <- p * (rad / d)
  return(p + ori)
}



# ===   7.4.2  Force directed layout: iterations

# Ready to run the FDL schedule

# In general, the number of iterations should be on the order of the number
# of nodes to consider
nIterate <- 500

stepSize <- 0.01  # move 1% of the force-vector magnitude in each iteration
kA <- 0.1         # attractive force constant

kRmin <- 0        # turn repulsive forces off initially
kRmax <- 0.004    #
kRwait <- 0.5     # wait fraction of steps long before ramping up kR
kRramp <- 0.2     # ramp up kR over this fraction of steps

Fmax <- 1.0       # in case of near overlap, repulsive forces can get
                  # VERY large. To be safe, we bound the magnitude of
                  # the force vector.

stepSums <- 0     # We don't really need this here - but keeping
                  # track of the movement during iterations allows you
                  # to see how fast it converges.


for (i in 1:nIterate) {  # iterate steps

  pBar(i, nIterate)


  # calculate current kR
  if (i < kRwait * nIterate) {  # kR is minimum
    kR <- kRmin
  } else if (i > (kRwait + kRramp) * nIterate) {   # kR is maximum
    kR <- kRmax
  } else {  # kR is being ramped up
    kR <- kRmin + (((i - (kRwait*nIterate))/(kRramp*nIterate)) * (kRmax-kRmin))
  }

  # prepare index of movable points in random order
  iMovable <- sample(which(fdl$isMovable))

  # one iteration: update all movable points
  for (idx in iMovable) {

    p <- c(fdl$x[idx], fdl$y[idx])   # fetch point coordinates ...
    pName <- fdl$name[idx]           # ... and name

    # collect the points that exert a force on p. These are:
    #   the gene it is linked to (same Name) ...
    sel <- which(fdl$name == pName & !(fdl$isMovable))
    vPoints <- matrix(c(fdl$x[sel], fdl$y[sel]), nrow = 1)

    #   ... and the protein(s) it makes functional interactions with (same
    #       row in myFIedges).
    FInames <- c(myFIedges$b[myFIedges$a == pName],
                 myFIedges$a[myFIedges$b == pName])
    sel <- which(fdl$name %in% FInames & fdl$isMovable) # movables only
    vPoints <- rbind(vPoints, matrix(c(fdl$x[sel], fdl$y[sel]),
                                     nrow = length(sel)))

    # calculate the force vectors
    thisForce <- fA(p, vPoints, kA) + fR(p, vPoints, kR)

    # bound the force if necessary
    magF <- sqrt(thisForce[1]^2 + thisForce[1]^2)
    if (magF > Fmax) {
      thisForce <- thisForce * (Fmax / magF)
    }

    # update the position
    p <- p + (stepSize * thisForce)
    stepSums[i] <- stepSums[i] + (stepSize * magF)

    # project it back on the circle
    p <- toCircle(p, CHR20ORI, CHR20RAD * PROTSCALE)

    # update the coordinates
    fdl[idx, c("x", "y")] <- p
  }
}

# plot(1:nIterate, log10(stepSums))



# ===   7.4.3  Force directed layout: coordinate update

# Done calculating the layout. Now write the coordinates back to the
# shapes:

# drop all rows of not movable points from fdl
fdl <- fdl[fdl$isMovable, ]

for (i in 1:length(myShapes)) {

  if (myShapes[[i]]$category == "FDLmovable") {
    thisName <- myShapes[[i]]$name                              # get name
    sel <- which(fdl$name == thisName)                          # find its row
    myShapes[[i]]$centre <- as.numeric(fdl[sel, c("x", "y")])   # update centre
  }
}

# Define Gene-Protein edges
# These edges will be drawn as straight lines, with the color of the
# GO code.

for (i in 1:nrow(myGPedges)) {
  gName <- myGPedges$g[i]
  pName <- myGPedges$p[i]
  for (j in 1:length(myShapes)) {
    if (myShapes[[j]]$class == "gene" && myShapes[[j]]$name == gName) {
      gXY <- myShapes[[j]]$centre
      next
    }
    if (myShapes[[j]]$class == "protein" && myShapes[[j]]$name == gName) {
      pXY <- myShapes[[j]]$centre
      next
    }
  }
  myShapes[[length(myShapes)+1]] <- list(type = "line",
                                         class = "connector",
                                         p1 = gXY,
                                         p2 = pXY,
                                         stroke=myGOcolours[myGPedges$type[i]],
                                         sw = myGPedges$weight[i])
}


# ==   7.5  Layout functional interaction edges  ===============================

# Define Functional interaction edges

# These edges will be drawn as cubic bezier curves, with control points on a
# line to the circle centre.


FIRAD <- CHR20RAD * PROTSCALE * 0.5  # radius of circle on which we place the
                                     # bezier control points
FICOL <- "#377ACD3F"  # 75% transparent light blue
FIWEIGHT <- 1.5       # Stroke weight fo functional interactions (in points)

for (i in 1:nrow(myFIedges)) {
  aName <- myFIedges$a[i]
  bName <- myFIedges$b[i]
  for (j in 1:length(myShapes)) {
    if (myShapes[[j]]$class == "protein" && myShapes[[j]]$name == aName) {
      aXY <- myShapes[[j]]$centre
    }
    if (myShapes[[j]]$class == "protein" && myShapes[[j]]$name == bName) {
      bXY <- myShapes[[j]]$centre
    }
  }

  cpaXY <- toCircle(aXY, CHR20ORI, FIRAD)
  cpbXY <- toCircle(bXY, CHR20ORI, FIRAD)

  myShapes[[length(myShapes)+1]] <- list(type = "cubic",
                                         class = "connector",
                                         p1 = aXY,
                                         p2 = cpaXY,
                                         p3 = cpbXY,
                                         p4 = bXY,
                                         stroke = FICOL,
                                         sw = myGPedges$weight[i])
}



# Done. All shapes are defined.


# =    8  PLOT  ================================================================
# cf. https://www.w3.org/TR/SVG

# ==   8.1  Compute scale and translation  =====================================

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



# ==   8.2  Write SVG header  ==================================================
mySVG <- SVGheader()
mySVG <- c(mySVG, SVGdefinePage(Xpx, Ypx))


# ==   8.3  Render all elements  ===============================================
#
for (i in 1:length(myShapes)) {

  mySVG <- c(mySVG, SVGrenderElement(myShapes[[i]],
                                     sc = sXY,
                                     tr = tXY,
                                     Y = Ypx))
}


# ==   8.4  Write SVG footer  ==================================================
mySVG <- c(mySVG, SVGfooter())


# =    9  FINISH  ==============================================================

# Write the SVG to file
writeLines(mySVG, con = SVGFILE)

# Open the SVG in the default browser to visualize
system(sprintf("open -a \"Google Chrome\" %s", SVGFILE))   # For MacOS
# Windows ???
# Linux ???



# ====  TESTS  =================================================================
# ...





# [END]
