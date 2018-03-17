import svgwrite
from setUpCh import chSetUp

CHROMLEN = 64444167 #nts on chromosome 20
START = 50          #left margin offset
DISPLEN = 400       #lengh to draw the chromosome
AXIS = 50           #top margin offset for chromosome
TITLEAXIS = 30      #top margin offset for title (auto centered)
BARHEIGHT = 20      #size of gene rects

#gene location based print - linear
#draw the chromosome as a line and add rectangles to represent gene occupancy

def linearGeneDraw(chDF, fName):

    #initialize svg drawing
    dwg = svgwrite.Drawing(filename=fName)

    #add title
    dwg.add(dwg.text("CHR 20", text_anchor="middle", insert=(START + DISPLEN / 2, TITLEAXIS)))

    #add chromosome line... want to scale everything by the display length
    dispScale = DISPLEN / CHROMLEN

    #draw the chromosome line according to left and top margins
    chEnd = CHROMLEN * dispScale + START

    #draw the chromosome line
    ch = dwg.add(dwg.g(id="ch", stroke="black", stroke_width=3))
    ch.add(dwg.line(start=(START, AXIS),
                    end=(chEnd, AXIS)))

    #add rects for each gene
    #insert upper left corner 1/2 the rect height above the chromosome line
    genes = dwg.add(dwg.g(id="genes", stroke="black", stroke_width=0.01))

    for g in chDF.index:
        gStart = chDF.loc[g].start * dispScale + START  #chose where to draw the rect
        gEnd = chDF.loc[g].end * dispScale + START      #chose where to end the rect

        #draw rect
        genes.add(dwg.rect(insert=(gStart, (AXIS-(BARHEIGHT/2))),
                           size=((gEnd-gStart), BARHEIGHT),
                           fill=chDF.loc[g].colour))
    #save the svg
    dwg.save()


if __name__ == '__main__':
    print("setting up...")
    ch20 = chSetUp()
    print("set up complete\nrendering svg...")
    linearGeneDraw(ch20, "linDraw.svg")
