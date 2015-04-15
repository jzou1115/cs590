import argparse
import random

def getAtoms(pdb):
    atoms=[]
    for line in pdb:
        tokens=line.split()
        if(tokens[0].strip() == "ATOM"):
            atoms.append(line)
    return atoms


#def checkAtoms(atoms):
#    i = int(atoms[0].split()[1].strip())
#    for a in atoms:
#        print a
#        if int(a.split()[1].strip()) != i:
#            return False
#        else:
#            i+=1
#    return True
 
def checkAtoms(atoms):
    return True

def newPDB(pdb, f):
    outfile = open(f, "w")
    for l in pdb:
        tokens=line.split()
        if tokens[0].strip() != "TER":
            outfile.write(l)
    outfile.close() 

def writeSysCfg(name, start, end, mut):
    outfile = open("System.cfg", "w")
    outfile.write("pdbname "+name+"\n")
    outfile.write("strand0 "+str(start)+" "+str(end)+"\n")
    outfile.write("strandMutNums "+mut+"\n")
    mutations = getMuts(start, end, int(mut))
    outfile.write("strandMut0 "+mutations+"\n")
    outfile.write("strandAA0  true\n")
    outfile.write("strandRotTrans false\n")
    outfile.write("onlySingleStrand 0\n")
    outfile.write("numOfStrands 1")
    outfile.close()

def getMuts(start, end, mut):
    mutations=[]
    while len(mutations)< int(mut):
        i = random.randrange(start, end)
        if i not in mutations:
            mutations.append(i)
    ret=""
    for m in mutations:
        ret = ret+ str(m)+" "
    return ret.strip()

def writeDeeCfg(name, mut):
    outfile = open("DEE.cfg", "w")
    outfile.write("runName "+name+"\n")
    outfile.write("numMaxMut "+mut+"\n")
    outfile.write("doMinimize false\n")
    outfile.write("minEnergyMatrixName min.dat")
    outfile.close()

def main():
    parser= argparse.ArgumentParser(description="Get file")
    parser.add_argument("pdb")
    parser.add_argument("mut")
    args= parser.parse_args()

    pdb= open(args.pdb, "r").readlines()
    atoms= getAtoms(pdb)
    if checkAtoms(atoms):
	#write new pdf file w/0 TER
	# newPDB(pdb, args.pdb+"terRemoved.pdb")

        #write system.cfg
        start = int(atoms[0].split()[5])
        end = int(atoms[len(atoms)-1].split()[5])
        writeSysCfg(args.pdb, start, end, args.mut)

        #write dee.cfg
        writeDeeCfg(args.pdb, args.mut)
    

main()
