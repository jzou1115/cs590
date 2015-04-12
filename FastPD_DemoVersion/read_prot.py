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

#lcosts=array('l', [int(round(costs[i][i])) for i in xrange(len(costs[0]))])
lcosts=array('l', [0 for i in xrange(numpoints*numlabels)])
print("lcosts: "+str(len(lcosts)))

for i in xrange(len(costs)):
	costs[i][i]=0
pairs=array('l', [item for sublist in [[(i,j) for j in xrange(i+1, len(costs[0]))] for i in xrange(len(costs[0]))][:-1] for pair in sublist for tup in (pair, pair[::-1]) for item in tup])
print("pairs: "+str(len(pairs)))

<<<<<<< HEAD
num=input("how long to wait: ")

labels=label_file.readline().split('\t')
costs=[map(float,line.split('\t')[:-1]) for line in cost_file]

numpoints=bytes(len(costs[0]))
numlabels=bytes(len(labels))
lcosts=bytes(array('d', [costs[i][i] for i in xrange(len(costs[0]))]))
pairs=bytes(array('d', [item for sublist in [[(i,j) for j in xrange(i+1, len(costs[0]))] for i in xrange(len(costs[0]))][:-1] for pair in sublist for item in pair]))
numpairs=bytes(len(pairs))
dist=bytes(array('d', [dist for lst in costs for dist in lst]))
max_iters=bytes((30))
wcosts=bytes(array('d', [1 for i in xrange(len(pairs))]))

data = [numpoints, numlabels, lcosts, numpairs, pairs, dist, max_iters, wcosts]
for item in data:
	out_file.write(item)
out_file.close()


=======
numpairs=len(pairs)/2
print("numpairs: "+str(numpairs))

dist = [0 for i in xrange(numlabels*numlabels)]
for i in xrange(len(costs)):
	for j in xrange(len(costs[i])):
		dist[i*numlabels + j] = int(round(costs[i][j]))
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
out_file.write(struct.pack('=l', numlabels))
out_file.write(struct.pack('=l', numpairs))
out_file.write(struct.pack('=l', max_iters))
out_file.write(struct.pack('='+'l'*len(lcosts), *lcosts))
out_file.write(struct.pack('='+'l'*len(pairs), *pairs))
out_file.write(struct.pack('='+'l'*len(dist), *dist))
out_file.write(struct.pack('='+'l'*len(wcosts), *wcosts))
out_file.close()


r=open('/home/aditya/git/cs590/FastPD_DemoVersion/matrices.bin', 'rb')
print(struct.unpack('l',r.read(8)))
print(struct.unpack('l',r.read(8)))
print(struct.unpack('l',r.read(8)))
print(struct.unpack('l',r.read(8)))
r.close()
>>>>>>> origin/adi
