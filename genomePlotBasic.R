# genomePlotBasic.R
#
# Purpose:  Demo a genome plot of genes on a line.
#
#
# Version:  1.0
# Date:     2018 03 18
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Dependencies:
#           readr package
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Version history:
#   1.0  Final version for Biohacks 2018
#   0.3  Added an introduction section and more comments
#   0.2  Updated gene data
#   0.1  Derived from genomePlotDemo.R
#
# ToDo:
#    - ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                  Line
#TOC> ------------------------------------------------------
#TOC>   1        INTRODUCTION                             47
#TOC>   2        PARAMETERS                              100
#TOC>   3        PACKAGES AND FUNCTIONS                  125
#TOC>   4        PROCESS                                 134
#TOC>   4.1        READ SOURCE DATA                      139
#TOC>   4.2        INITIALIZE DATA STRUCTURES            147
#TOC>   4.3        ANNOTATE                              165
#TOC>   4.4        LAYOUT                                172
#TOC>   4.5        PLOT                                  286
#TOC>   4.5.1          Compute scale and translation     289
#TOC>   4.5.2          Write SVG header                  316
#TOC>   4.5.3          Render all elements               321
#TOC>   4.5.4          Write SVG footer                  332
#TOC>   5        FINISH                                  336
#TOC>
#TOC> ==========================================================================


# =    1  INTRODUCTION  ========================================================

# This script demonstrates a data-driven visualization of genome scale data. It
# is a very basic example, but written in a modular way that should make it easy
# to extend into something more interesting. Extensions of the code are provided
# in genomePlotDemo.R and genomePlotIntermediate.R (There is no "advanced"
# version - that would be your code).

# The code proceeds through five steps:

#  1 -  Read the source data:
#         A single data file containing gene annotations is read into a
#         data frame. Extensions would read any additional kind of data
#         that contributes to the final image at this point.
#
#  2 -  Initialize data structures:
#         In this step, information objects are defined. In this demo we define
#         only gene information, and we store it in a data frame, one row per
#         gene. Extensions would add non-gene information like general
#         phentotypes, categories of ethnicity, geaographical distribution,
#         a map of subcellular locations, tissue diagrams, explanatory text
#         and legends, network data - whatever the visualization should include.
#
#  3 -  Annotate:
#         In this demo, no further annotation is performed. Extensions would
#         correlate, compile and contrast data, to populate attributes of our
#         data structures with values.
#
#  4 -  Layout:
#         In this step data is mapped to visuals. We call these "shapes", and
#         we build them from the items in our data structures. In this
#         script we define a single line to represent Chr 20, add a text label,
#         and plot each gene as a coloured rectangle on the line, where
#         position represents the chromosomal location. This step implements
#         the main vision about how the data should be visualized and we
#         expect that it will absorb the majority of your effort. Extensions
#         would include arranging shapes in a variety of different basic
#         layouts, adding lines and curves that connect shapes to symbolize
#         processes and relationships, implementing masks, and layers,
#         optimizing placement of shapes with algorithms such as
#         force-directed layout, or hierarchical trees,
#         transforming sections of the layout with local deformations such as
#         perspective and hyperbolic lenses and much, much more.
#
#  5 -  Plot:
#         In this step the shapes that have been defined previously are
#         rendered to SVG markup that can be saved and displayed. A basic
#         set of functions have been provided for this and we expect that
#         any extensions you need are easy to add. Speak to the mentors for
#         advice.




# =    2  PARAMETERS  ==========================================================

# Code should not contain "magic numbers". Constants that we need are
# defined and commented here.

CHR20LENGTH   <- 64444167  # basepairs of the chromosome we are working with

DATAFILE <- "Chr20GeneData.tsv"  # Chromosome data input. See README-DATA and
                                 # prepareGenomeData.R for a description of the
                                 # contents, and the code that produces it from
                                 # database sources.
DATAFILE <- "Chr20Reduced.tsv" # alternative file with reduced number of lines

NACOLOUR <- "#AAAAAA"            # Neutral grey for NA attributes

SVGFILE <- "chromosome.svg"  # Filename for the output we produce

# UTPoster prints from 24" x 36" all the way to 60" x 300".
# Let's assume legal paper size for this demo, landscape orientation, and
# subtract a 1" margin on both sides.
# PAGEWIDTH  <- ( 14.0 - 2) * 2.54    # in cm
# PAGEHEIGHT <- (  8.5 - 2) * 2.54    # in cm

PAGEWIDTH  <- (  8.5 - 2) * 2.54 
PAGEHEIGHT <- ( 8.0 - 2) * 2.54
RESOLUTION <- 150                  # pixels per 2.54 cm



# =    3  PACKAGES AND FUNCTIONS  ==============================================
#
# All required packages and functions are loaded from the source file below.
# You can inspect/copy/modify the source code there.
setwd("/home/dennis/BioHacks/2018-Challenge/")
source("genomePlotFunctions.R")



# =    4  PROCESS  =============================================================

# This demo code will plot genes as rectangles on a line and color the boxes.


# ==   4.1  READ SOURCE DATA  ==================================================

# read_tsv() is from the readr package. It is similar to base R's read.delim()
# function, but more modern.

myData <- read_tsv(DATAFILE)


# ==   4.2  INITIALIZE DATA STRUCTURES  ========================================

# There are many possibilities to store the data for the objects we will analyze
# and draw. Here we take a very simple approach and store gene-level data in one
# data frame.

# Entity annotations: data for each gene
myGenes <- data.frame(sym = myData$sym,            # Gene symbols
                      start = myData$start+0.05*CHR20LENGTH,        # start
                      end = myData$end+0.05*CHR20LENGTH,            # end
                      strand = myData$strand,      # strand
                      GOid = myData$GO_P,          # GO annotation for "Process"
                      stringsAsFactors = FALSE)



# ==   4.3  ANNOTATE  ==========================================================


# (No computed annotations in this basic demo.)



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
# chromosome backbone as a line and we will add the chromosome name as a text
# element.

# At first, we are only concerned with relative positions and we will layout
# shapes into an arbitrary canvas. Later we will map this into page coordinates.
# The chromosome line will be placed on the X axis, start at (0, 0) and run
# to (1, 0)

# Chromosome backbone:

# CHR20FROM <- c(0.05, 0.0)
# CHR20TO   <- c(1.05, 0.0)
CHR20FROM <- c(0.5, 0.05)
CHR20TO   <- c(0.5, 1.05)

myShapes[[1]] <- list(type = "line",
                      p1 = CHR20FROM,
                      p2 = CHR20TO,
                      stroke = "#4499AA",  # colour of line
                      sw = 7.0)            # stroke-width

# Next we add some descriptive text:

myShapes[[1]] <- list(type = "text",
                      text = "",
                      centre = c(0.1, 0.1),
                      size = 48,           # points
                      font = "Times",
                      fill = "#33AAFF")

# Next we add the genes. We will draw them as rectangles. genes on the (+)
# strand will go above the circle, genes on the (-) strand will go below the
# circle. We will place them on the circle at their fractional position on the
# chromosome. We will also give them colour:

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

  # The centre of the rectangle is placed above or below the line, depending
  # on the strand. Calculate height first:
  # thisHeight <- 0.02 * (CHR20TO[1] - CHR20FROM[1])   # 1% of chromosome length
  thisWidth <- 0.02 * (CHR20TO[2] - CHR20FROM[2])

  thisHeight <- abs(myGenes$start[i] - myGenes$end[i]) / CHR20LENGTH
  #thisHeight <- abs(myGenes$start[i] - myGenes$end[i]) / CHR20LENGTH
  
  #thisCentre <- c(mean(c(myGenes$start[i], myGenes$end[i])) / CHR20LENGTH,  # x
  #             0.0 + (myGenes$strand[i] * thisHeight * 0.5) )                # y
  thisCentre <- c(0.5 + (myGenes$strand[i] * thisWidth * 0.5),  # x
                mean(c(myGenes$start[i], myGenes$end[i])) / CHR20LENGTH )           # y

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
                                           centre = thisCentre,
                                           w = thisWidth,
                                           h = thisHeight,
                                           ang = 0,
                                           fill = myFill,
                                           stroke = myFill,
                                           sw = 0.5)  # points
}

# add the final strand again, since we want that on top
myShapes[[length(myShapes) + 1]] <- list(type = "line",
                      p1 = CHR20FROM,
                      p2 = CHR20TO,
                      stroke = "#4499AA",  # colour of line
                      sw = 7.0)            # stroke-width
# Done. All shapes are defined.


# ==   4.5  PLOT  ==============================================================
# cf. https://www.w3.org/TR/SVG

# ===   4.5.1  Compute scale and translation

# Caution: the SVG coordinate system has its origin (0, 0) in the TOP LEFT
# corner, positive X goes right, and positive Y goes down. Here we define the
# necessary scaling and translation.


# range of coordinates:
#dX <- (CHR20TO[1] - CHR20FROM[1]) * 1.1
dX <- (CHR20TO[2] - CHR20FROM[2]) * 1.1
# Given the range that needs to fit on the page, we can compute the scale:

sXY <- RESOLUTION * (PAGEWIDTH / 2.54) / dX

# We compute the dimensions of the page in pixels ...

Xpx <- RESOLUTION * (PAGEWIDTH  / 2.54)
Ypx <- RESOLUTION * (PAGEHEIGHT / 2.54)

# And we compute a translation: for this demo, we move the chromosome
# to the middle of the page, and 1% to the right:

tXY <- c(0.01 * (CHR20TO[2] - CHR20FROM[2]),
         Ypx / 2)



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
