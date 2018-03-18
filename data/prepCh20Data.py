import os
import pandas as p

# Read file obtained via custom download from ensembl biomart
# Custom download:
# Filters: Chromosome/scaffold: 20
#          With HGNC Symbol ID(s): Only
# Attrib.: Gene stable ID
#          Gene start (bp)
#          Gene end (bp)
#          Strand
#          GOSlim GOA Accession(s)
#          GOSlim GOA Description
#          HGNC symbol

def prepData():
    #load the biomart data
    rawFile = os.path.abspath("mart_export.txt")
    rawData = p.read_table(rawFile, header=0)

    #don't need the gene ID column but we do want a chromosome column. We can mutate this in place.
    rawData["Gene stable ID"] = 20

    #rename columns
    rawData = rawData.rename(columns={
                            "Gene stable ID": "chr",
                            "Gene start (bp)": "start",
                            "Gene end (bp)": "end",
                            "Strand": "strand",
                            "GOSlim GOA Accession(s)": "GOAid",
                            "GOSlim GOA Description": "GOAdescr",
                            "HGNC symbol": "sym"})

    #set sym to be the index and remove the sym col
    rawData = rawData.set_index(rawData.sym)
    rawData = rawData.drop("sym", axis=1)

    #make GOA frequency table
    GOAfreq = rawData.GOAid.value_counts()

    #remove duplicate rows
    rawData = rawData.drop_duplicates()

    #choose the most specific (least frequent) GOA
    for symbol in rawData.index.unique():
        #if there are multiple, select row with most specific GOA
        if p.Series(rawData.loc[symbol].GOAid).count() > 1:
            topGOA = GOAfreq[rawData.loc[symbol].GOAid].idxmin()
            #drop all other rows for that symbol -> keep all those that are not the symbol of interest OR have the topGOA
            rawData = rawData[((rawData.index != symbol) | (rawData.GOAid == topGOA))]

    #sort output by HUGO symbol
    rawData = rawData.sort_index()

    #write to tsv
    rawData.to_csv("chr20_data.tsv", sep="\t")

if __name__ == '__main__':
    prepData()