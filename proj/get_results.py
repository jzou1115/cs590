import sys

COST_FILE = open(sys.argv[1], 'r')
RES_FILE = open(sys.argv[2], 'r')
OUT_FILE = open(sys.argv[3], 'w')

COUNTS = [int(x) for x in COST_FILE.readline().split('\t')]
LABELS = COST_FILE.readline().split('\t')

for val in RES_FILE:
    OUT_FILE.write(LABELS[int(val)]+'\n')

COST_FILE.close()
RES_FILE.close()
OUT_FILE.close()
