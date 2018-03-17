#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 11:35:34 2018

@author: dennis
"""

import pandas as pd
import numpy as np
import math
import os
import json
    
# parameters and data
filename = "Chr20GWAStraits.tsv"
output = "edges.json"
sep = "\t"

# manually defined dictionary of classification mappings
class_dict = {"asthma":"physical",
              "arthritis":"physical",
              "diabetes":"physical",
              "narcolepsy":"physical",
              "cancer":"physical",
              "leukemia":"physical",
              
              "depression": "mental",
              "anger": "mental",
              "schizophrenia": "mental",
              "bipolar": "mental",
              "attention": "mental",
              "loneliness": "mental",
              
              "head": "head",
              "face": "head",
              "facial": "head",
              "brain": "head",
              "intelligence": "head",
              "well-being": "head",
              "memory": "head",
              "hair": "head",
              
              "heart": "heart",
              
              "blood": "blood",
              "artery": "blood",
              
              "lung": "lung",
              "pulmonary": "lung",
              
              "breast": "breast",
              
              "waist": "waist",
              "obesity": "waist",
              "bmi": "waist",
              "body mass": "waist",
              
              "liver": "liver",
              "kidney": "liver",
              
              "height": "height",
              
              }

if __name__ == "__main__":
    os.chdir(".")
#    interactions = pd.read_csv("Chr20FuncIntx.tsv", sep="\t")
    
#    traits = pd.read_csv("Chr20GWAStraits.tsv", sep="\t")

    
    # Manually read in file. Split at specified separator, clean any line breaks, and force lowercase
    trait_list = []
    with open(filename, "r") as f:
        data = f.readlines()
        for line in data:
            temp = line.split(sep)
            temp[1] = temp[1].lower().strip("\n")
            trait_list.append(temp)

    # skip header
    trait_list = trait_list[1:]
    
    # create empty sets and dictionaries
    # We're using sets, since it automatically cancels the problem of duplicate edges.
    mapping = set()
    score = {}
    
    # go over all the examples we have and find whether keywords appear
    for key in class_dict.keys():
        for j, el in enumerate(trait_list):
            if key in el[1]:
                tup = (el[0], class_dict[key])
                
                # add an intensity score, which is linearly increased
                if (tup in mapping):
                    score[tup] += 1
                else:
                    score[tup] = 1
                mapping.add(tup)
                
    # create output file

    # combine to list, so we can put it to JSON
    temp_list = []
    for el in mapping:
        # line layout: source (Gene), target (bodypart), score
        temp_list.append([el[0], el[1], score[el]])
        
    json_string = json.dumps(temp_list)
    
    with open(output, "w") as f:
        f.write(json_string)
    


