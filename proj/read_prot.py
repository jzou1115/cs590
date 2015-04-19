"""
read_prot_long
"""

import struct
import sys

LABEL_FILE = sys.argv[1]
OUT_FILE = sys.argv[2]

COST_FILE = open(LABEL_FILE, 'r')
OUT_FILE = open(OUT_FILE, 'wb')

#read in rotamer COUNTS per residue
COUNTS = [int(x) for x in COST_FILE.readline().split('\t')]
#get label-rotamer correspondence
LABELS = COST_FILE.readline().split('\t')
#unaltered pairwise energy matrix
COSTS = [map(float, line.split('\t')) for line in COST_FILE]

NUMPOINTS = len(COUNTS)
NUMLABELS = len(LABELS)

#get label COSTS from energy matrix
INTERNALS = [COSTS[i][i] for i in xrange(len(COSTS))]
LABELINGS = [INTERNALS for i in xrange(NUMPOINTS)]
for i in xrange(len(LABELINGS)):
    before = sum(COUNTS[:i])
    count = COUNTS[i]
    l = [2147483647 for x in xrange(before)]
    l.extend(INTERNALS[before:before+count])
    l.extend([2147483647 for x in xrange(before+count, NUMLABELS)])
    LABELINGS[i] = l


#move label COSTS to a 1-D array
LCOSTS = [0 for i in xrange(NUMLABELS*NUMPOINTS)]
for i in xrange(len(LABELINGS)):
    for k in xrange(len(LABELINGS[i])):
        LCOSTS[k*NUMPOINTS+i] = LABELINGS[i][k]

#this is the most ridiculous loop comprehension I have ever written
PAIRS = [item
    for sublist in [[(i, j)
        for j in xrange(i+1, NUMPOINTS)]
        for i in xrange(len(COSTS[0]))][:-1]
    for pair in sublist
    for item in pair]
NUMPAIRS = len(PAIRS)/2

#remove label cost from pairwise COSTS
for i in xrange(len(COSTS)):
    COSTS[i][i] = 0

#move label COSTS to a 1-D array
DIST = [0 for i in xrange(NUMLABELS*NUMLABELS)]
for i in xrange(len(COSTS)):
    for j in xrange(len(COSTS[i])):
        DIST[j*NUMLABELS + i] = COSTS[i][j]

WCOSTS = [1 for i in xrange(NUMPAIRS)]
MAX_ITERS = 50

PACK_LCOSTS = '=' + 'd'*len(LCOSTS)
PACK_PAIRS = '=' + 'l'*len(PAIRS)
PACK_DIST = '=' + 'd'*len(DIST)
PACK_WCOSTS = '=' + 'l'*len(WCOSTS)

OUT_FILE.write(struct.pack('=l', NUMPOINTS))
OUT_FILE.write(struct.pack('=l', NUMPAIRS))
OUT_FILE.write(struct.pack('=l', NUMLABELS))
OUT_FILE.write(struct.pack('=l', MAX_ITERS))
OUT_FILE.write(struct.pack(PACK_LCOSTS, *LCOSTS))
OUT_FILE.write(struct.pack(PACK_PAIRS, *PAIRS))
OUT_FILE.write(struct.pack(PACK_DIST, *DIST))
OUT_FILE.write(struct.pack(PACK_WCOSTS, *WCOSTS))

OUT_FILE.close()
COST_FILE.close()
