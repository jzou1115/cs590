

with open('../../OSPREY_v2.2beta_Feb2015/dataFiles/LovellRotamer.dat.vol', 'r') as f:
    for line in f:
        print len(line.split(' '))-1
