#!/usr/bin/env python

import pandas as pd
import os

## -----------------------------------------------------------------------------------------
# Read in files
## -----------------------------------------------------------------------------------------

# import paths
depthpath = os.environ.get("DEPTH")
depth = pd.read_csv(depthpath, sep='\t', lineterminator='\n', names = ["chr","position","depth"])

## -----------------------------------------------------------------------------------------
## Step 1: Look at the coverage before and after sites 
## -----------------------------------------------------------------------------------------
n = depth.shape[0]

def condition(x): 
  return x > 2

if (n == 4):
  depth["site"] = ['A', 'B', 'C', 'D']
  siteTest1 = sum(condition(x) for x in depth.depth[depth['site'] == "A"])
  siteTest2 = sum(condition(x) for x in depth.depth[depth['site'] == "B"])
  siteTest3 = sum(condition(x) for x in depth.depth[depth['site'] == "C"])
  siteTest4 = sum(condition(x) for x in depth.depth[depth['site'] == "D"])
  outside = siteTest1 + siteTest4
  inside = siteTest2 + siteTest3
elif (n == 3):
  depth["site"] = ['A', 'B', 'C']
  siteTest1 = sum(condition(x) for x in depth.depth[depth['site'] == "A"])
  siteTest2 = sum(condition(x) for x in depth.depth[depth['site'] == "B"])
  siteTest3 = sum(condition(x) for x in depth.depth[depth['site'] == "C"])
  outside = siteTest1 + siteTest3
  inside = siteTest2
else: 
  print("error")

## -----------------------------------------------------------------------------------------
## Step 5: Check integration type
## -----------------------------------------------------------------------------------------

if ((inside > 0) and (outside > 0)):
    print("linear_duplication")
elif ((inside > 0) and (outside == 0)):
    print("ecDNA_integration")
elif ((inside == 0) and (outside > 0)):
    print("deletion_integration")
else:
    print("undertermined_double_integration")
