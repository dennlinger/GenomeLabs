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


if __name__ == "__main__":
    os.chdir(".")
    interactions = pd.read_csv("Chr20FuncIntx.tsv", sep="\t")
#    pass


