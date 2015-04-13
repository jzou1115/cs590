from array import array
import struct
import numpy
import sys

out_file=open('/home/aditya/git/cs590/FastPD_DemoVersion/matrices.bin', 'wb+')

numpoints=2
numlabels=5
labelings=[[1, 2, 9001, 9001, 9001], [9001, 9001, 1, 2, 3]]

lcosts=[0 for i in xrange(numlabels*numpoints)]
for i in xrange(len(labelings)):
	for k in xrange(len(labelings[i])):
		lcosts[k*numpoints+i]=labelings[i][k]
numpairs=1
pairs=[0, 1]

costs=[[0, 9001, 3, 2, 1], [9001, 0, 3, 2, 1], [2, 1, 0, 9001, 9001], [2, 1, 9001, 0, 9001], [2, 1, 9001, 9001, 0]]
dist=[0 for i in xrange(25)]
for i in xrange(len(costs)):
	for j in xrange(len(costs[i])):
		dist[j*numlabels + i]=costs[i][j]

wcosts=[1]
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
