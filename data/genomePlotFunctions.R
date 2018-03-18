# genomePlotFunctions.R
#
# Functions for genome plot starter code.
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author:  where not otherwise stated:
#          Boris Steipe <boris.steipe@utoronto.ca>
# Notes:
#
# ==============================================================================
#


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                            Line
#TOC> ----------------------------------------------------------------
#TOC>   1        PACKAGES                                           40
#TOC>   2        FUNCTIONS TO INITIALIZE DATA STRUCTURES            48
#TOC>   3        ANNOTATION FUNCTIONS                               51
#TOC>   4        LAYOUT FUNCTIONS                                   66
#TOC>   4.01       ang2rad()                                       110
#TOC>   4.02       coord2circle()                                  128
#TOC>   5        PLOTTING FUNCTIONS                                156
#TOC>   5.01       SVGheader()                                     158
#TOC>   5.02       SVGdefinePage()                                 176
#TOC>   5.03       SVGdrawCircle()                                 200
#TOC>   5.04       SVGdrawText()                                   231
#TOC>   5.05       SVGdrawLine()                                   269
#TOC>   5.06       SVGdrawCubic()                                  296
#TOC>   5.07       SVGdrawRect()                                   327
#TOC>   5.08       st()                                            375
#TOC>   5.09       SVGrenderElement()                              401
#TOC>   5.10       SVGfooter()                                     482
#TOC> 
#TOC> ==========================================================================


# =    1  PACKAGES  ============================================================

if (!require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}


# =    2  FUNCTIONS TO INITIALIZE DATA STRUCTURES  =============================


# =    3  ANNOTATION FUNCTIONS  ================================================

linMap <- function(x, low, high) {
  # map a vector of numbers x to values between min and max. NA values
  # are ignored.

  d <- high - low
  x <- x  - min(x, na.rm = TRUE)
  x <- x * (d / max(x, na.rm = TRUE))
  x <- x + low
  return(x)
}



# =    4  LAYOUT FUNCTIONS  ====================================================


category2colour <- function(Cs,
                            col,
                            B = 1,            # bias
                            I = "spline",     # interpolate
                            A = FALSE) {      # use alpha channel too
  # Description:
  #    Define a named vector of colours for the categories in Cs according to
  #    the colour palette in pal.
  #
  # Arguments:
  #     Cs:          a vector items that can be used as rownames
  #     pal:  char   a palette of colour values
  #     B:    num    see bias in colorRampPalette() help
  #     I:    char   <"linear", "spline">
  #     A:    bool   return alpha channel values or not
  #
  # Details:
  #     Cs can have any atomic type. If pal has only one element, all the
  #     colours are the same. Two elements will create a gradient, more
  #     then two elements will make a colour spectrum. The other parameters
  #     have reasonable defaults but can be adjusted.
  #
  # Value:
  #     A named character vector of colour values, where the names are
  #     the elements of Cs.

  fCol <- colorRampPalette(colors = col,
                           bias = B,
                           interpolate = I,
                           alpha = A)

  cV <- fCol(length(Cs))
  names(cV) <- Cs

  return(cV)
}





# ==   4.01  ang2rad()  ========================================================
ang2rad <- function(a) {
  # Description:
  #    Convert a rotation angle a in degrees from vertical in clockwise
  #    direction to radians.
  #
  # Arguments:
  #     a:   num   rotation angle in degrees from 0 at (0, 1), clockwise.
  #
  # Value:
  #     An angle in radians.

  x <- ((2 * a) / 360) * pi  # degrees to radians
  x <- -x                    # change direction of rotation
  x <- x + (pi / 2)          # add 90°
  return(x)
}

# ==   4.02  coord2circle()  ===================================================
coord2circle <- function(coord, l, ori, rad) {
  # Description:
  #     Convert linear coordinates on a line of length l to
  #     positions on a circle centred on ori with radius rad. Return the
  #     coordinates and the rotation angle. We define 0° to be at the top, and
  #     the positive direction is clockwise.
  #
  # Arguments:
  #     coord:   num   chromosome coordinates
  #     l:       num   chromosome length
  #     ori:     num   (x, y) coordinates of the circle centre
  #     rad:     num   circle radius
  #
  # Value:
  #     A three element vector of (x, y) coordinates on the circle and the
  #     rotation angle in degree, where 0 is a vertical line and the
  #     direction of rotation is clockwise.

  rot <- (coord / l) * 360
  th  <- ang2rad(rot)
  x <- (cos(th) * rad) + ori[1]
  y <- (sin(th) * rad) + ori[2]

  return(c(x, y, rot))
}


# =    5  PLOTTING FUNCTIONS  ==================================================

# ==   5.01  SVGheader()  ======================================================
SVGheader <- function() {
  #
  # Description:
  #     Write a SVG header.
  #
  # Arguments:
  #     None
  #
  # Value:
  #     String. SVG header

  s <-    "<?xml version=\"1.0\" standalone=\"no\"?>"
  s[2] <- paste("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">")
  return(s)
}

# ==   5.02  SVGdefinePage()  ==================================================
SVGdefinePage <- function(w, h) {
  #
  # Description:
  #     Define pahe size and viewbox
  #
  # Arguments:
  #     w:   Page width (pixels)
  #     h:   Page height (pixels)
  #
  # Value:
  #     String. SVG page definitions.

  s <- paste("<svg",
             sprintf("width=\"%f\"",  w),
             sprintf("height=\"%f\"", h),
             sprintf("viewbox=\"0 0 %f %f\"", w, h),
             "version=\"1.1\"",
             "xmlns=\"http://www.w3.org/2000/svg\">",
             collapse = " ")
  return(s)
}


# ==   5.03  SVGdrawCircle()  ==================================================
SVGdrawCircle <- function(cx, cy, r,
                          fill = "#FFFFFF",
                          stroke= "#000000",
                          sw = 1.0) {
  #
  # Description:
  #     Write SVG code to draw a circle.
  #
  # Arguments:
  #     cx      num   x coordinate of centre
  #     cy      num   y coordinate of centre
  #     r       num   radius
  #     fill    char  fill colour, default white
  #     stroke  char  outline colour, default black
  #     sw      char  stroke width, default 1.0 point
  #
  # Value:
  #     String. SVG commands.

  s <- paste( sprintf("<circle cx=\"%f\"", cx),
              sprintf("cy=\"%f\"", cy),
              sprintf("r=\"%f\"", r),
              sprintf("fill=\"%s\"", fill),
              sprintf("stroke=\"%s\"", stroke),
              sprintf("stroke-width=\"%f\" />", sw),
              collapse = " ")
  return(s)
}


# ==   5.04  SVGdrawText()  ====================================================
SVGdrawText <- function(x, y,
                        font,
                        size,
                        fill = "#000000",
                        text) {
  #
  # Description:
  #     Write SVG code to place text.
  #
  # Arguments:
  #     x       num   x coordinate of centre
  #     y       num   y coordinate of centre
  #     font    char  a valid font-family
  #     size    num   font-size in points
  #     fill    char  fill colour, default black
  #     text    char  the text string
  #
  # Details:
  #     The text is centred on (x, y).
  #
  # Value:
  #     String. SVG commands.

  s <- paste( sprintf("<text x=\"%f\"", x),
              sprintf("y=\"%f\"", y),
              "alignment-baseline=\"middle\"", # this centers the string on
              "text-anchor=\"middle\"",        # the (x, y) coordinates
              sprintf("font-family=\"%s\"", font),
              sprintf("font-size=\"%f\"", size),
              sprintf("fill=\"%s\">", fill),
              collapse = " ")
  s <- c(s, text)
  s <- c(s, "</text>")
  return(s)
}


# ==   5.05  SVGdrawLine()  ====================================================
SVGdrawLine <- function(x1, y1, x2, y2,
                        stroke = "#000000",
                        sw = 1.0) {
  #
  # Description:
  #     Write SVG code to draw a line.
  #
  # Arguments:
  #     x1, y1, x2, y2  num   start and end coordinates of the line
  #     stroke          char  line colour, default black
  #     sw              num   stroke width, default 1.0 point
  #
  # Value:
  #     String. SVG commands.

  s <- paste( sprintf("<line x1=\"%f\"", x1),
              sprintf("y1=\"%f\"", y1),
              sprintf("x2=\"%f\"", x2),
              sprintf("y2=\"%f\"", y2),
              sprintf("style=\"stroke:%s;", stroke),
              sprintf("stroke-width:%f\" />", sw),
              collapse = " ")
  return(s)
}


# ==   5.06  SVGdrawCubic()  ===================================================
SVGdrawCubic <- function(x1, y1, x2, y2, x3, y3, x4, y4,
                         fill = "none",
                         stroke = "#000000",
                         sw = 1.0) {
  #
  # Description:
  #     Write SVG code to draw a cubic Bezier curve.
  #
  # Arguments:
  #     x1, y1          num   start coordinates of the curve
  #     x2, y2, x3, y3  num   control points
  #     x4, y4          num   end coordinates of the line
  #     stroke          char  line colour, default black
  #     sw              num   stroke width, default 1.0 point
  #
  # Value:
  #     String. SVG commands.

  s <- paste( sprintf("<path d=\""),
              sprintf("M %f %f ", x1, y1),
              sprintf("C %f %f %f %f %f %f ", x2, y2, x3, y3, x4, y4),
              sprintf("Z\" "),
              sprintf("fill=\"%s\" ", fill),
              sprintf("stroke=\"%s\" ", stroke),
              sprintf("stroke-width=\"%f\" />", sw),
              collapse = " ")
  return(s)
}


# ==   5.07  SVGdrawRect()  ====================================================
SVGdrawRect <- function(x, y, w, h, ang,
                        fill = "#FFFFFF",
                        stroke = "none",
                        sw = 0.0) {
  #
  # Description:
  #     Write SVG code to draw a rectangle.
  #
  # Arguments:
  #     x, y    num   coordinates of centre
  #     w, h    num   rectangle width and height
  #     ang     num   rotation angle of the rectangle
  #     fill    char  fill colour, default white
  #     stroke  char  outline colour, default: no outline
  #     sw      num   stroke width, default 1.0 point
  #
  # Details:
  #     The rectangle is contained in a <g></g> element that defines the
  #     rotation by ang degrees, clockwise, as a transformation.
  #
  # Value:
  #     String. SVG commands.

  crX <- x + (w / 2)
  crY <- y + (h / 2)
  #transform=\"",
  #sprintf("rotate(%f, %f, %f)", ang, crX, crY),\
  s <- paste( "<g>",
              collapse = " ")
  s <- c(s, paste( sprintf("<rect x=\"%f\"", x),
                   sprintf("y=\"%f\"", y),
                   sprintf("width=\"%f\"", w),
                   sprintf("height=\"%f\"", h),
                   sprintf("fill=\"%s\"", fill),
                   sprintf("stroke=\"%s\"", stroke),
                   ifelse(stroke != "none",
                          sprintf("stroke-width=\"%f\"", sw),
                          ""),
                   "/>",
                   collapse = " "))
  s <- c(s, "</g>")

  return(s)
}


# ==   5.08  st()  =============================================================
st <- function(xy, s, t, yPage) {
  # Description:
  #   Simple scaling and translation
  #
  # Arguments:
  #   xy:    num  two-element vector of the coordinates to be scaled
  #   s:     num  scale
  #   t:     num  two element tranlation vector
  #   yPage  num  height of the page in pixel
  #
  # Details:
  #   Note that the SVG coordinate system has its origin top-left.
  #   This functioncould be extended to rotations,
  #   perspective skews, projections etc. etc.
  #
  # Value:
  #     transformed coordinate pair.

  xy <- s * xy                      # scale
  xy <- c(xy[1], yPage - xy[2])     # flip
  xy <- xy + c(t[1], - t[2])        # translate
  return(xy)
}


# ==   5.09  SVGrenderElement()  ===============================================

SVGrenderElement <- function(li, sc, tr, Y) {
  #
  # Description:
  #     Render a shape to an SVG plot command
  #
  # Arguments:
  #     li:   list element that contains the shape type and attributes
  #     sc:   num   scale
  #     tr:   num   translation vector
  #     Y:    num   page height (in pixel)
  #
  # Details:
  #     This sends the shapes to the appropriate SVG function with a minimum
  #     of preprocessing. See actual function calls for details.
  #
  # Value:
  #     String. SVG commands.

  if (li$type == "circle") {
    xy <- st(li$centre, sc, tr, Y)
    s <- SVGdrawCircle(cx = xy[1],
                       cy = xy[2],
                       r = li$radius * sc,
                       fill = li$fill,
                       stroke = li$stroke,
                       sw = li$sw)
  } else if (li$type == "text") {
    xy <- st(li$centre, sc, tr, Y)
    s <- SVGdrawText(x = xy[1],
                     y = xy[2],
                     font = li$font,
                     size = li$size,
                     fill = li$fill,
                     text = li$text)
  }  else if (li$type == "line") {
    xy1 <- st(li$p1, sc, tr, Y)
    xy2 <- st(li$p2, sc, tr, Y)
    s <- SVGdrawLine(x1 = xy1[1],
                     y1 = xy1[2],
                     x2 = xy2[1],
                     y2 = xy2[2],
                     stroke = li$stroke,
                     sw = li$sw)
  }  else if (li$type == "cubic") {
      xy1 <- st(li$p1, sc, tr, Y)
      xy2 <- st(li$p2, sc, tr, Y)
      xy3 <- st(li$p3, sc, tr, Y)
      xy4 <- st(li$p4, sc, tr, Y)
      s <- SVGdrawCubic(x1 = xy1[1],
                        y1 = xy1[2],
                        x2 = xy2[1],
                        y2 = xy2[2],
                        x3 = xy3[1],
                        y3 = xy3[2],
                        x4 = xy4[1],
                        y4 = xy4[2],
                        stroke = li$stroke,
                        sw = li$sw)
    }  else if (li$type == "rect") {
      xCorner <- li$centre[1] - (li$w / 2)
      yCorner <- li$centre[2] + (li$h / 2)
      xy <- st(c(xCorner, yCorner), sc, tr, Y)
      s <- SVGdrawRect(x = xy[1],
                       y = xy[2],
                       w = li$w * sc,
                       h = li$h * sc,
                       ang = li$ang,
                       fill = li$fill,
                       stroke = li$stroke,
                       sw = li$sw)
    } else {
    stop(sprintf("Unsupported SVG element type \"%s\"", li$type))
  }

  return(s)
}



# ==   5.10  SVGfooter()  ======================================================
SVGfooter <- function() {
  #
  # Description:
  #     Return closing tag for svg file.
  #
  # Arguments:
  #     None
  #
  # Value:
  #     The string "</svg>"

  return("</svg>")
}


# [END]

