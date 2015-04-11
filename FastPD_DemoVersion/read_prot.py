from array import array
import struct
import sys

# variable assignments are as follows:
#	 - numpoints: number of MRF nodes

# 	- numlabels: number of labels

#	 - lcosts: 1-dimensional array of size numlabels*numpoints containing the 
#	   label costs (i.e., the MRF singleton-potentials). The label cost
#	   for the k-th label at the i-th node is assumed to be given by
#	   lcosts[k*numpoints+i].

#	 - numpairs: number of MRF edges

#	 - int *pairs: 1-dimensional array of size 2*numpairs containing the 
#	   nodes' indices for each MRF edge. The indices for the two nodes of 
#	   the i-th edge are assumed to be given by: pairs[2*i], pairs[2*i+1].
#	   (Note that the nodes' indices start from 0).

#	 - Real *dist: the distance function used for defining the MRF pairwise
#	   potentials. This is an 1-dimensional array of size numlabels*numlabels.
#	   The distance D(x_p,x_q) is assumed to be given by: dist[x_q*numlabels+x_p].
#	   (Label indices start from 0).
  
#	 - max_iters: maximum number of outer iterations of the FastPD algorithm

#	 - wcosts: 1-dimensional array of size numpairs containing the weights w_{pq}
#	   used in the MRF pairwise potentials. wcosts[i] is the weight corresponding 
#	   to the i-th MRF edge.

in_file=open('/home/aditya/Desktop/proteins/fastpd/osprey.txt', 'r')
out_file=open('/home/aditya/Desktop/proteins/fastpd/matrices.bin', 'w+b')

numpoints = int(f.readline()).to_bytes(4, sys.biteorder)
numlabels = int(f.readline()).to_bytes(4, sys.biteorder)
lcosts = array('d', [float(item) for item in f.readline().split(',')])
numpairs = int(f.readline()).to_bytes(4, sys.biteorder)
pairs = array('d', [int(item) for item in f.readline().split(',')])
dist = array('d', [float(item) for item in f.readline().split(',')])
max_iters = (30).to_bytes(4, sys.bite_order)
wcosts = array('d', [1 for i in xrange(numpairs)])

data = [numpoints, numlabels, lcosts, numpairs, pairs, dist, max_iters, wcosts]
for item in data:
	item.to_file(out_file)
out_file.close()