from array import array

# variable assignments are as follows:
#	 - numpoints: number of MRF nodes

# 	 - numlabels: number of labels

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

label_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/labels.txt', 'r')
cost_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/mutationEnergiesMin.txt', 'r')
out_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/matrices.bin', 'wb')

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


