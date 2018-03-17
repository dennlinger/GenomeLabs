import os
import pandas as p

def chSetUp():
    #load in data - see prepCh20Data.py
    ch20file = os.path.abspath("chr20_data.tsv")
    ch20 = p.read_table(ch20file, header=0)

    #make rows selectable based on HUGO symbol
    ch20 = ch20.set_index("sym")

    #assign different colours to each unique GOAid
    GOAs = p.DataFrame()
    GOAs["GOAid"] = p.Series(ch20.GOAid.unique())

    #get colours
    colfile = os.path.abspath("svgCols.txt")
    cols = p.read_table(colfile, header=None)

    #assign colours and place into ch20 dataframe
    GOAs["colour"] = cols[0:len(GOAs.GOAid)]
    ch20 = ch20.join(GOAs.set_index("GOAid"), on="GOAid")

    return ch20