#!/usr/bin/env python

import re

# DEFINE FUNCTIONS

def shared(x, y):
	sh = 0
	grX = x.split("|")
	grY = y.split("|")
	
	for i in range(len(grX)):
		grX[i] = set(grX[i].split(","))
	for j in range(len(grY)):
		grY[j] = set(grY[j].split(","))
		
		
	for xXx in grX:
		for yYy in grY:
			if xXx.issubset(yYy):
				sh += 1
				break
				
	return sh

def scoreConflict(indA, patA, patterns):
    vrai = True

    dividand = 0

    pg = patA.split("|")
    score = 0.0
    pats2 = patterns[:]
    overall_rank = 0.0
    conflict_score = 0.0
    print "PATS2:"
    print pats2
    print "del pats2[%d]" % indA
    del pats2[indA]
    print pats2
    for bb, patB in enumerate(pats2):
        if re.search("\|", patB):  # don't include constant sites!!
            dividand += 1
            compB = patB.split("|")
            scB = 0.0
            for group in compB:
                cont = 1
                if vrai:
                    tax = group.split(",")
                    ch1 = tax[0]
                    for gp in pg:
                        if ch1 in gp:
                            for ch2 in tax[1:]:
                                if ch2 not in gp and cont:
                                    scB = scB + 1.0
                                    cont = 0
            print "%s v %s = %f" % ( patA, patB, (scB/len(compB)) )
            score += (scB/len(compB))

    conflict_score = score/dividand
    print "conflict_score: %s vs %s\t%f/%d\t%f" % (patA, patB, score, dividand, conflict_score)
	
    print "******************************"
    return conflict_score
    
def uniqify(seq):
    seen = {}
    result = []
    for item in seq:
        if item in seen : continue
        seen[item] = 1
        result.append(item)
    return result


def getPattern(site, unknown):
    considered = []
    pattern = []
    for x in range(len(site)):
        if site[x] not in unknown:
            if site[x] in considered:
                pattern[considered.index(site[x])].append(str(x))
	    else:
                considered.append(site[x])
		pattern.append([str(x)])


    patStr = "|".join([",".join(g) for g in pattern])

    return patStr


def jumblePattern(site, unknown):
	import random

	siteJ = ""
	while site:
		pos = random.randrange(len(site))
		siteJ += site[pos]
		site = site[:pos] + site[(pos+1):]

	return getPattern(siteJ, unknown)


def DNAdetect(seq):
    seq = seq.upper()
    oLen = float(len(seq))
    seq_C = ""

    seq_C = seq.replace("A", "")
    seq_C = seq_C.replace("C", "")
    seq_C = seq_C.replace("G", "")
    seq_C = seq_C.replace("T", "")
    
    nLen = float(len(seq_C))
    perc = (nLen/oLen)*100

    if perc < 20.0:
        return "DNA"
    else:
        seq_C = seq.replace("0", "")
        seq_C = seq_C.replace("1", "")
        if len(seq_C) == 0:
            return "standard"
        else:
            return "protein"
        

def histogram(num_list, name_list):
    upper = float(max(num_list))
    pad = len(str(upper))
    parts = []
    for i in range(1,61):
        parts.append((upper/60)*i)

    for m, n in enumerate(num_list):
        pr = name_list[m] + "|"
        low = 0.0
        if n == 0:
            pr = pr + " "*61
        for p, hi  in enumerate(parts):
            if n > low and n <= hi:
                pr = pr + "="*(p+1) + (" "*(60 - p))
                break
            low = hi
        print "[" + pr + "|" + str(n) + " "*(pad-len(str(n))) + "]"


def FastaParse(file_name):
    file = open(file_name)
    names = []
    seqs = []
    for line in file:
	    if ">" in line:
		    names.append(line[1:].strip())
		    seqs.append("")
	    else:
		    seqs[-1] += line.strip()

    ret = [names, seqs]
    return ret


def printHelp():
    print """
****************
TIGER Help:
****************

TIGER: Tree-Independent Generation of Evolutionary Rates

(Developed by Carla Cummins in the lab of James Mc Inerney, NUI Maynooth, Co. Kildare, Ireland)

-Options:

-in        Specify input file. File must be in FastA format and must be aligned prior.
           Datasets with uneven sequence lengths will return an error.

-v         Returns current TIGER version.

-f         Changes output formatting options.
           -f s: sorts sites depending on their agreement score
           -f r: displays rank values rather than bin numbers
           -f s,r: displays sorted ranks (*Be sure to put only a "," NO SPACE!)
           
           Default prints bin numbers unsorted.

-b         Set the number of bins to be used.
           -b <int>: Sites will be placed into <int> number of bins. <int> is a whole number.

           Default is 10

-rl       A list of the rate at each site may be optionally written to a specified
          file. 
          -rl <file.txt> : writes list of the rates at each site to file.txt.

-ptp      Specifies that a PTP test should be run. *Note: this option has a huge 
          effect on running time!

-z        Number of randomisations to be used for the PTP test. 
          -z <int>: each site will be randomised <int> times. <int> is a whole number.

	 Default is 100

-p       Specify p-value which denotes significance in PTP test.
         -p <float>: site will be denoted as significant if p-value is better than <float>.
                     <float> is a floating point number.

         Default is 0.05

-pl      Write a list of p-values to a specified file.
         -pl <file.txt>: writes list of p-values for each site to file.txt.

-u       Specify unknown characters in the alignment. Unknown characters are omitted from 
         site patterns and so are not considered in the analysis.
         -u ?,-,*: defines ?, - and * as unknown characters. (*Be sure to put only a comma
                   between characters, NO SPACE!!)
         
         Default is ? only


-Examples:
     
1.   ./TIGER -in ExampleFile.aln -f s,r -v -rl rate_list.txt

     This will run the software on "ExampleFile.aln", with sorted ranks included in the output.
     The variability measure for each site will be displayed and a list of the rates at (unsorted)
     sites will be written to the file "rate_list.txt".
  
2.   ./TIGER -in ExampleFile.aln -ptp -r 1000 -p 0.01 -u ?,*

     This will run the software on the file "ExampleFile.aln" with a PTP test. Sites will be
     randomised 1,000 times and pass the test if their p-value is <0.01. All ? and * characters
     encountered in the alignment will be ommitted from the analysis.
   
     """
