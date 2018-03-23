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


## Further Ideas
We were confronted with several setbacks during the development process, but have some more ideas which could (some quite easily) be realized for further development:

* Extension to Gene-specific traits (Human Protein Atlas (HPA) description)
* Extension to biomedically motivated regions (adding keywords,...)
* Displaying more information, potentially from an API-accessible database
* View that represents the chromosome more accurately (e.g., circular view around person, or highlighted based on current cursor position (zooming view))
* Give indicator whether those are beneficial or malevolent traits (i.e. cancer vs longer lifetime, ...; indicated by some HPA values)

