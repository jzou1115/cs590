from array import array
import struct
import numpy
import sys

label_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/labels.txt', 'r')
cost_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/mutationEnergiesMin.txt', 'r')
out_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/matrices.bin', 'wb')

labels=label_file.readline().split('\t')
costs=[map(float,line.split('\t')[:-1]) for line in cost_file]

###############################################################################################################################

#numpoints=2
#numlabels=5
#labelings=[[1, 2, 9001, 9001, 9001], [9001, 9001, 1, 2, 3]]

#lcosts=[0 for i in xrange(numlabels*numpoints)]
#for i in xrange(len(labelings)):
#        for k in xrange(len(labelings[i])):
#                lcosts[k*numpoints+i]=labelings[i][k]
#numpairs=1
#pairs=[0, 1]

#dist=[0 for i in xrange(25)]
#for i in xrange(len(costs)):
#        for j in xrange(len(costs[i])):
#                dist[j*numlabels + i]=costs[i][j]

#wcosts=[1]
#max_iters=30

################################################################################################################################
numpoints=2
numlabels=len(labels)

internals=[costs[i][i] for i in xrange(len(costs))]
label1 = internals[0:3]
label1.extend([0 for x in xrange(3, len(internals))])
label2 = [0 for x in xrange(3)]
label2.extend(internals[3:])
labelings=[label1, label2]

lcosts=[0 for i in xrange(numlabels*numpoints)]
for i in xrange(len(labelings)):
	for k in xrange(len(labelings[i])):
		lcosts[k*numpoints+i]=labelings[i][k]

pairs=[item for sublist in [[(i,j) for j in xrange(i+1, numpoints)] for i in xrange(len(costs[0]))][:-1] for pair in sublist for item in pair]
numpairs=len(pairs)/2

for i in xrange(len(costs)):
	costs[i][i] = 0

dist=[0 for i in xrange(numlabels*numlabels)]
for i in xrange(len(costs)):
	for j in xrange(len(costs[i])):
		dist[j*numlabels + i] = costs[i][j]

wcosts=[1 for i in xrange(numpairs)]
max_iters=30


################################################################################################################################



#numpoints=2
#numlabels=len(labels)

#maxint=2147483647
#internals=[int(round(costs[i][i])) for i in xrange(len(costs[0]))]
#lcosts=[item for item in internals*len(costs[0])]

#for i in xrange(len(costs)):
#	costs[i][i]=0
#pairs=array('l', [item for sublist in [[(i,j) for j in xrange(i+1, len(costs[0]))] for i in xrange(len(costs[0]))][:-1] for pair in sublist for item in pair])


#numpairs=len(pairs)/2

#dist = [0 for i in xrange(numlabels*numlabels)]
#for i in xrange(len(costs)):
#	for j in xrange(len(costs[i])):
#		dist[i*numlabels + j] = int(round(costs[i][j]))
#		if (j==0 and i!=0):
#			dist[i*numlabels + j] = 2147483647

#max_iters=30

#wcosts=array('l', [1 for i in xrange(numpairs)])

# The input binary file is assumed to contain all the necessary data for
# defining the energy of a discrete MRF. 
# More specifically, it should contain the following variables according 
# to the order indicated below:
#  *  numpoints
#  *  numlabels
#  *  numpairs
#  *  max_iters
#  *  lcosts
#  *  pairs
#  *  dist
#  *  wcosts

out_file.write(struct.pack('=l', numpoints))
out_file.write(struct.pack('=l', numpairs))
out_file.write(struct.pack('=l', numlabels))
out_file.write(struct.pack('=l', max_iters))
out_file.write(struct.pack('='+'l'*len(lcosts), *lcosts))
out_file.write(struct.pack('='+'l'*len(pairs), *pairs))
out_file.write(struct.pack('='+'l'*len(dist), *dist))
out_file.write(struct.pack('='+'l'*len(wcosts), *wcosts))
out_file.close()
