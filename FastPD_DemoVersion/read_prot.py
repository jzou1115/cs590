import struct

label_file=open('labels.txt', 'r')
cost_file=open('mutationEnergiesMin.txt', 'r')
out_file=open('matrices.bin', 'wb')

#read in rotamer counts per residue
counts=map(int, label_file.readline().split('\t'))
#get label-rotamer correspondence
labels=label_file.readline().split('\t')
#unaltered pairwise energy matrix
costs=[map(float,line.split('\t')[:-1]) for line in cost_file]

numpoints=len(counts)
numlabels=len(labels)

#get label costs from energy matrix
internals=[costs[i][i] for i in xrange(len(costs))]
labelings=[internals for i in xrange(numpoints)]
for i in xrange(len(labelings)):
	before=sum(counts[:i])
	count=counts[i]
	l = [9001 for x in xrange(before)]
	l.extend(internals[before:before+count])
	l.extend([9001 for x in xrange(before+count, numlabels)])
	labelings[i]=l

#move label costs to a 1-D array
lcosts=[0 for i in xrange(numlabels*numpoints)]
for i in xrange(len(labelings)):
	for k in xrange(len(labelings[i])):
		lcosts[k*numpoints+i]=labelings[i][k]

#this is the most ridiculous loop comprehension I have ever written
pairs=[item for sublist in [[(i,j) for j in xrange(i+1, numpoints)] for i in xrange(len(costs[0]))][:-1] for pair in sublist for item in pair]
numpairs=len(pairs)/2

#remove label cost from pairwise costs
for i in xrange(len(costs)):
	costs[i][i] = 0

#move label costs to a 1-D array
dist=[0 for i in xrange(numlabels*numlabels)]
for i in xrange(len(costs)):
	for j in xrange(len(costs[i])):
		dist[j*numlabels + i] = costs[i][j]

wcosts=[1 for i in xrange(numpairs)]
max_iters=30

out_file.write(struct.pack('=l', numpoints))
out_file.write(struct.pack('=l', numpairs))
out_file.write(struct.pack('=l', numlabels))
out_file.write(struct.pack('=l', max_iters))
out_file.write(struct.pack('='+'l'*len(lcosts), *lcosts))
out_file.write(struct.pack('='+'l'*len(pairs), *pairs))
out_file.write(struct.pack('='+'l'*len(dist), *dist))
out_file.write(struct.pack('='+'l'*len(wcosts), *wcosts))

out_file.close()
label_file.close()
cost_file.close()
