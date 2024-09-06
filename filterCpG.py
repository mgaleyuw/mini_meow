#!/usr/bin/env python

#takes three positional arguments: first argument is region, second argument
#is input file. Outputs cpgIslands.argv1.CGonly.bed in directory
#supplied by arvgv3.

import numpy as np
import sys

np.set_printoptions(threshold=100000)

region =sys.argv[1]
infile=sys.argv[2]
outname = sys.argv[3]

try:
    f = open(infile, "r")
except:
    raise Exception("Input not found")
flines = f.read()
f.close()

frecords = flines.split(">")[1:]
g = open(outname, "w")
#h = open(tempname2, "w")
for fx in frecords:
    fxstart = int(fx.split("\n")[0].split(":")[-1].split("-")[0])
    fxtext = fx.split("\n")[1].upper()
    fxtextCG = np.array(list(fxtext.replace("CG", "XX")))
    fxname=region
    if region=="NONE":
        fxname=fx.split("\n")[0].split(":")[0]
    pos = np.arange(fxstart, fxstart+len(fxtext))
    cpgPos = pos[(fxtextCG=="X")]
    labels=np.array([fxname]*len(cpgPos))
    g.write(np.array2string(np.array([labels, cpgPos]).transpose(), separator="\t").replace(" ", "").replace("[", "").replace("]", "").replace("'", ""))
    g.write("\n")

g.close()
