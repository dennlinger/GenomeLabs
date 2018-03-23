# GenomeLabs

## Human Gene Atlas
This was the winning entry of the BCB BioHacks 2018 hackathon, hosted at the University of Toronto.<br/>

We developed a simple visualization of matchings of GWAS traits of genes represented in Chromosome 20, and linked them to specific parts of the body.

## Implementation
The main implementation is the visualization, which is delivered in a D3 sankey environment. The preprocessing of the data was done in Python, based on the provided files.

### Preprocessing
To generate the necessary bipartite graph network, we created a manually curated list of keywords, for which we looked up the trait description according to the Chr20GWAStraits.tsv file.
Any time a given key word would appear in the traits description, we added a link between the respective gene and keyword. We additionally (linearly) increase a weight of how many times a certain keyword and gene would appear together, which allows us to model the graph dynamically based on "connection strength". <br/>

Furthermore, we restricted our visualization to Genes that would appear in the list of extracted connections, which reduces the number of genes to around 170 of the total of 500 genes in the file.

We then manually parsed the relevant information into JSON.
That includes the edges, which we generated, as well as the gene information available form the Chr20GeneData.tsv file. Due to time constraints, we were not able to include the relevant GWAS data, which would have been a nice addition, to highlight where we found that information.

## Implementation
The JS-based implementation builds upon JQuery, and D3 - specifically sankey graphs. We built the basis on top of a tutorial we found online [http://bl.ocks.org/d3noob/5028304].


## Problems encountered
Some of the problems where the following:
* Opacity overwrite in sankey. Since we want to highlight links both during hovering, as well as hovering the corresponding nodes, we had a problem in overwriting CSS-defined opacity. Could maybe avoided by defining opacity through sankey directly
* Positioning/Ordering of the genes. Sankey is giving us a hard time defining the exact position/size of the graph.
* Adding relevant information to the file. JQuery is not made for appending a lot of information easily. Using node.js instead could simplify some of this.
* Responsive design for the canvas. The size and positioning make use of weirdly chosen divisions, which had to be manually adjusted. For now, the canvas has a fixed size.


## Further Ideas
We were confronted with several setbacks during the development process, but have some more ideas which could (some quite easily) be realized for further development:

* Extension to Gene-specific traits (Human Protein Atlas (HPA) description)
* Extension to biomedically motivated regions (adding keywords,...)
* Displaying more information, potentially from an API-accessible database
* View that represents the chromosome more accurately (e.g., circular view around person, or highlighted based on current cursor position (zooming view))
* Give indicator whether those are beneficial or malevolent traits (i.e. cancer vs longer lifetime, ...; indicated by some HPA values)


If time is given, it might be more helpful to move the idea to a more sustainable architecture, which can be extended in a more modular fashion; this would require the transition of some of the sankey code, since it does some aspects (Curvature, color choice) already pretty neatly.

