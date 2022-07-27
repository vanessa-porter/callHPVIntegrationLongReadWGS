import pandas as pd
import os

## -----------------------------------------------------------------------------------------
# Read in files
## -----------------------------------------------------------------------------------------

# import paths
path = os.environ.get("DIR")

## -----------------------------------------------------------------------------------------
# Read in test1 file
## -----------------------------------------------------------------------------------------

t1 = pd.read_csv((path + "/intTest1.txt"), sep='\t', lineterminator='\n', names = ["result"])

repeat = t1.result[0]
integ = t1.result[1]

# check for files
test2 = os.path.exists(path + "/intTest2.txt")
test3 = os.path.exists(path + "/intTest3.txt")

if (test2 == True) and (test3 == True):
    t2 = pd.read_csv(path + "/intTest2.txt", sep='\t', lineterminator='\n', names = ["result"])
    twoBreak = t2.result[0]
    t3 = pd.read_csv(path + "/intTest3.txt", sep='\t', lineterminator='\n', names = ["result"])
    ecDNA = t3.result[0]
    if (ecDNA != "ecDNA_conditions_not_met") and (twoBreak != "ecDNA_integration"):
        print("complex_ecDNA_integration")
    elif (twoBreak == "linear_duplication"):
        print("duplication_integration")
    elif (twoBreak == "deletion_integration"):
        print("deletion_integration")
    elif (twoBreak == "ecDNA_integration"):
        print("ecDNA_integration")
    else:
        print("undetermined_integration")

elif (test2 == False) and (test3 == True):
    t3 = pd.read_csv(path + "/intTest3.txt", sep='\t', lineterminator='\n', names = ["result"])
    ecDNA = t3.result[0]
    if (repeat == "centromere_repeats") or (repeat == "telomere_repeats"):
        print("repeat_integration")
    elif (ecDNA != "ecDNA_conditions_not_met"):
        print("complex_ecDNA_integration")
    elif (integ == "translocation_integration"):
        print("translocation_integration")
    elif (integ == "multi-breakpoint_integration"):
        print("multi-breakpoint_integration")
    elif (integ == "unmatched_integration"):
        print("unmatched_integration")
    else:
        print("undetermined_integration")

elif (test2 == True) and (test3 == False):
    t2 = pd.read_csv(path + "/intTest2.txt", sep='\t', lineterminator='\n', names = ["result"])
    twoBreak = t2.result[0]
    if (twoBreak == "linear_duplication"):
        print("duplication_integration")
    elif (twoBreak == "deletion_integration"):
        print("deletion_integration")
    elif (twoBreak == "ecDNA_integration"):
        print("ecDNA_integration")
    else:
        print("undetermined_integration")

elif (test2 == False) and (test3 == False):
    if (repeat == "centromere_repeats") or (repeat == "telomere_repeats"):
        print("repeat_integration")
    elif (integ == "translocation_integration"):
        print("translocation_integration")
    elif (integ == "multi-breakpoint_integration"):
        print("multi-breakpoint_integration")
    elif (integ == "unmatched_integration"):
        print("unmatched_integration")
    else:
        print("undetermined_integration")