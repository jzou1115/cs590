from array import array
import struct
import numpy

label_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/labels.txt', 'r')
cost_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/mutationEnergiesMin.txt', 'r')
out_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/matrices.bin', 'wb')

labels=label_file.readline().split('\t')
costs=[map(float,line.split('\t')[:-1]) for line in cost_file]



numpoints=len(costs[0])
print("numpoints: "+str(numpoints))

numlabels=len(labels)
print("numlabels: "+str(numlabels))

internals=[int(round(costs[i][i])) for i in xrange(len(costs[0]))]
internals[0]=10000
lcosts=[item for item in internals*len(costs[0])]
print("lcosts: "+str(len(lcosts)))

for i in xrange(len(costs)):
	costs[i][i]=0
pairs=array('l', [item for sublist in [[(i,j) for j in xrange(i+1, len(costs[0]))] for i in xrange(len(costs[0]))][:-1] for pair in sublist for item in pair])
print("pairs: "+str(len(pairs)))
print(pairs)


numpairs=len(pairs)/2
print("numpairs: "+str(numpairs))

dist = [0 for i in xrange(numlabels*numlabels)]
for i in xrange(len(costs)):
	for j in xrange(len(costs[i])):
		dist[i*numlabels + j] = int(round(costs[i][j]))
		if (j==0 and i!=0):
			dist[i*numlabels + j] = 10000
print("dist: "+str(len(dist)))

max_iters=30
print("max_iters: "+str(max_iters))

wcosts=array('l', [1 for i in xrange(numpairs)])

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
