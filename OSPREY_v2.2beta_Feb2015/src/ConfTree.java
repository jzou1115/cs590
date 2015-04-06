
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.PriorityBlockingQueue;

/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.2 beta
	Copyright (C) 2001-2014 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
 */

///////////////////////////////////////////////////////////////////////////////////////////////
//	ConfTree.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 *
 * @author mhall44
 */
public class ConfTree implements Serializable {
    //This is meant to be a general representation of a tree for conformational search
    //by A* or conceivably some other branch-and-bound type thing
    
    int numTreeLevels;
    
    PriorityBlockingQueue<ConfTreeNode> curExpansion = new PriorityBlockingQueue<>();
    
    int numExpanded = 0;
    int numFS = 0;
    
    boolean fixedOrder = true;//no dynamic ordering
    
    //Each tree level initially has numSimpleOptions[level] options (e.g., rotamers at a given res position)
    //but we can also create compound options, which will be indexed in levelOptions.get(level)
    //and are a list of the simple options making up the compound option
    //(so for example levelOptions.get(res).get(rot) should just be {rot}, and then later we can have options
    //like {rot1,rot2})
    //in ConfTreeNodes, conf indicates the compound option index at each level
    //which is just the simple option index for a simple option
    //then levelOptions.get(res).get(numSimpleOptions[res]) will be allowing everything
    int numSimpleOptions[];
    ArrayList<ArrayList<ArrayList<Integer>>> levelOptions;
    
    
    //stuff specific to rotamer expansion
    ReducedEnergyMatrix pairwiseMinEnergyMatrix;
    PairwiseEnergyMatrix pemUnreduced;
    
    
    boolean[][][][][][] splitFlags = null;
    boolean[][][][][][][][][] tripleFlags = null;
    //WERE GONNA TRY TO USE THESE IN FSCOREPBLP...NOT USED IF NULL
    
    CETMatrix cetm;
    
    
    PairwiseEnergyMatrix rigidMtx;//rigid-rot matrix can be optionally provided for UB calc
    //RCs should parallel pemUnreduced

    public ConfTree(ReducedEnergyMatrix pairwiseMinEnergyMatrix, PairwiseEnergyMatrix pemUnreduced) {
        this.pairwiseMinEnergyMatrix = pairwiseMinEnergyMatrix;
        this.pemUnreduced = pemUnreduced;
        
        //construct other stuff based on reduced energy matrix
        numTreeLevels = pairwiseMinEnergyMatrix.numRes;
        numSimpleOptions = new int[numTreeLevels];
        for(int res=0; res<numTreeLevels-1; res++)
            numSimpleOptions[res] = pairwiseMinEnergyMatrix.resOffsets[res+1]-pairwiseMinEnergyMatrix.resOffsets[res];
        
        numSimpleOptions[numTreeLevels-1] = pairwiseMinEnergyMatrix.numRotTotal 
                - pairwiseMinEnergyMatrix.resOffsets[numTreeLevels-1];
        
        levelOptions = new ArrayList<>();
        
        for(int res=0; res<numTreeLevels; res++){
            levelOptions.add(new ArrayList<ArrayList<Integer>>());
            for(int rot=0; rot<numSimpleOptions[res]; rot++){
                ArrayList<Integer> singleton = new ArrayList<>();
                singleton.add(rot);
                levelOptions.get(res).add(singleton);
            }
            
            //also add totally unassigned option
            ArrayList<Integer> all = new ArrayList<>();
            for(int rot=0; rot<numSimpleOptions[res]; rot++)
                all.add(rot);
            levelOptions.get(res).add(all);
        }
    }
    
    
    
        
    //compute the lowest-energy conformation remaining in the tree (i.e., in the union of the nodes's conf spaces)
    public int[] doAStar (boolean run1, int numMaxChanges, int nodesDefault[], boolean prunedNodes[],
			StrandRotamers strandRot[], String strandDefault[][], int numForRes[], int strandMut[][], boolean singleSeq, 
			int mutRes2Strand[], int mutRes2MutIndex[]){


		if (run1) {//if this is the first run of A*, then we need to set-up the empty expansion queue

                        //start by just putting in all undecided conformation
                        int conf[] = new int[numTreeLevels];
                        for(int level=0; level<numTreeLevels; level++)
                            conf[level] = numSimpleOptions[level];
                        
                        ConfTreeNode newNode = makeNode(conf);
                        curExpansion.put(newNode);
		}

		boolean done = false;		
		//start forming the conformation by expanding the lowest-valued node and updating the expansion queue
		/*
		 * While conf not fully defined
		 * 	For the current minimum node in the queue
		 * 		For all possible nodes at the next level
		 * 			If the steric is allowed
		 * 				Compute the f(n) scores; Add to the queue
		 * 		Delete the expanded node from the queue
		 */
                ConfTreeNode expNode = null;
                
		while (!done) {	
			
			expNode = curExpansion.poll();//take out min node from queue
			
			if (expNode==null){//the queue is empty
                            int emptyConf[] = new int[numTreeLevels];
                            Arrays.fill(emptyConf,-1);
                            return emptyConf;//so return a sequence of -1's to flag the stop of the search
			}

                        //expNode = updateFScore(expNode);//If lazily tightening F-score, do this now
                        
                        ConfSpaceSplit curSplit = pickSplit(expNode);
                        
                        if(curSplit==null)//no split available: conformation fully defined
                            done = true;
                        else {
                            ArrayList<ConfTreeNode> children = getChildren(expNode, curSplit);
                            
                            for(ConfTreeNode child : children){

                                if(child!=null){
                                    boolean addNode = true;

                                    //NOT YET DOING MAX NUMBER OF MUTATIONS FOR SPLITBYSPLACK
                                    /*
                                    if(!singleSeq){
                                        if(countMutants(child,nodesDefault)>numMaxChanges)
                                            addNode = false;
                                    }*/

                                    if( addNode )
                                        curExpansion.put(child);
                                }
                            }
                            
                            numExpanded++;
                        }
		}

                System.out.println("ConfTree: A* returning conformation; lower bound = "+expNode.LB+" nodes expanded: "+numExpanded+" FS terms evaluated: "+numFS);
		for(int a : expNode.conf)
                    System.out.print(a+" ");
                System.out.println();
                
                //if we get here expNode.conf should be all simple options
                return expNode.conf;
    }
    
    
    
    int getFirstSplittableLevel(int[] conf){
        for(int level=0; level<numTreeLevels; level++){
            if(conf[level]>=numSimpleOptions[level])//level not fully assigned
                return level;
        }
        return numTreeLevels;//indicates no splittable levels
    }
    
    ConfSpaceSplit pickSplit(ConfTreeNode expNode){
        
        int firstSplittableLevel = getFirstSplittableLevel(expNode.conf);
        if(firstSplittableLevel==numTreeLevels)//all levels fully assigned
            return null;//no more splits possible
                
        
        if(fixedOrder){
            //want full split of firstSplittableLevel into simple options
            ArrayList<Integer> opt = new ArrayList<>();
            for(int simpleOpt=0; simpleOpt<numSimpleOptions[firstSplittableLevel]; simpleOpt++)
                opt.add(simpleOpt);
            return new ConfSpaceSplit(firstSplittableLevel,opt,this);
        }
        
        //if not fixed order let's use the split with best harmonic average of quick f score increases relative to start
        ArrayList<ConfSpaceSplit> possibleSplits = listPossibleSplits(expNode);
        
        double bestRise = 0;
        ConfSpaceSplit bestSplit = possibleSplits.get(0);
        double origScore = quickFScore(expNode.conf);
        
        for(ConfSpaceSplit split : possibleSplits){
            ArrayList<ConfTreeNode> children = getChildren(expNode,split);//children created with quick f score
            
            double recipSum = 0;
            int numChildren = 0;
            for(ConfTreeNode child : children){
                double nodeRise = child.LB - origScore;
                if(nodeRise>0){
                    recipSum += 1./nodeRise;
                    numChildren++;
                }
                else{//no rise: don't consider this split
                    //(unless all splits have no rise--then take split # 0 by default)
                    numChildren = 0;
                    break;
                }
            }
            
            double harmRise = 1./(recipSum/numChildren);
            
            if(harmRise>bestRise){
                bestRise = harmRise;
                bestSplit = split;
            }
        }
        
        return bestSplit;
    }
    
    
    ArrayList<ConfSpaceSplit> listPossibleSplits(ConfTreeNode expNode){
        //what splits can be made to the conformational space of this ConfTreeNode?
        //if we get here we can assume not fixed order
        
        ArrayList<ConfSpaceSplit> ans = new ArrayList<>();
        
        //controls
        boolean useSeqSplits = false;
        boolean binarySplitOffs = false;//split off one rot at a time
        
        for(int res=0; res<numTreeLevels; res++){
            int curOption = expNode.conf[res];
            
            
            /*if(curOption==-1&&useSeqSplits){//must split into different AA types
                ans.add( new ConfSpaceSplit(res,makeAATypeOptions(res)) );
                make AA type sets cacheable;
                also go on and try other splits if this ends up singleton;
            }
            else */if(curOption>=numSimpleOptions[res]){//not fully specified at this residue
                
                if(binarySplitOffs){//split off one simple option at a time
                    for(int simpleOption : levelOptions.get(res).get(curOption)){
                        
                        ArrayList<Integer> allButOne = new ArrayList<>();
                        for(int opt : levelOptions.get(res).get(curOption)){
                            if(opt!=simpleOption)
                                allButOne.add(opt);
                        }
                        int otherOption = makeCompoundOption(res,allButOne);
                        
                        ArrayList<Integer> binSplit = new ArrayList<>();
                        binSplit.add(simpleOption);
                        binSplit.add(otherOption);
                        
                        ans.add(new ConfSpaceSplit(res,binSplit,this));
                    }
                }
                else {//split fully into simple options
                    ArrayList<Integer> opt = new ArrayList<>();
                    for(int simpleOpt : levelOptions.get(res).get(curOption))
                        opt.add(simpleOpt);
                    ans.add(new ConfSpaceSplit(res,opt,this));
                }
                
            }
        }
        
        return ans;
    }
    
    
    
    int makeCompoundOption(int level, ArrayList<Integer> simpleOptionList){
        //we want the index for a compound option at the specified residue consisting of the
        //listed simple options.  Create it if needed
        //we assume no repeats in simpleOptionList
        //This can probably be sped up by hashing or something if necessary
        
        //first, see if there is already a compound option like this
        for(int opt=0; opt<levelOptions.get(level).size(); opt++){
            ArrayList<Integer> curOpt = levelOptions.get(level).get(opt);
            if(curOpt.size()==simpleOptionList.size()){
                boolean match = true;
                for(int simpleOpt : simpleOptionList){
                    if(!curOpt.contains(simpleOpt)){
                        match = false;
                        break;
                    }
                }
                if(match)
                    return opt;
            }
        }
        
        //compound option doesn't exist...create it
        levelOptions.get(level).add(simpleOptionList);
        return levelOptions.get(level).size()-1;
    }
    
    
    ConfTreeNode updateFScore(ConfTreeNode expNode){//If lazily tightening F-score, do this now
        
        throw new RuntimeException("NOT SUPPORTED YET...");
        
        /*if(useFitSeries){
            //We will only expand terms that have the FS term included
            //until our lowest lower bound contains an FS term, we will
            //compute this term for the lowest-bounded nodes we get
            //and then throw them back in the queue (this is to minimize the number of
            //times we have to compute the FS term)
            while(!expNode.FSTermIncluded){
                expNode.fScore += FSTerm(expNode);
                expNode.FSTermIncluded = true;
                curExpansion.put(expNode);
                expNode = curExpansion.poll();
            }
            //Now expNode has an FS term so we can expand it
        }
        
        return expNode;*/
    }
    

    double fScoreRigid(int conf[], PairwiseEnergyMatrix pem){
        double f = 0;
            
        for(int level=0; level<numTreeLevels; level++){

            double bestE = Double.POSITIVE_INFINITY;//best energy (w.r.t. rotamer choices)
            //for the current level, including interactions with all
            //subsequent residues

            for( int rot : levelOptions.get(level).get(conf[level]) ){
                int index1 = pairwiseMinEnergyMatrix.resOffsets[level] + rot;
                bestE = Math.min(bestE, rotE(level,index1,conf,pem) );
            }                

            f += bestE;
        }
        
        return f;
    }
    
    
    double rotE(int level, int index1, int[] conf, PairwiseEnergyMatrix pem){
        //contribution to quickFScore for given rotamer (reduced, long index) at given level
        //this version of rotE is to use an alternate matrix!!!!
        int aa1 = pairwiseMinEnergyMatrix.indicesEMatrixAA[index1];
        int rc1 = pairwiseMinEnergyMatrix.indicesEMatrixRot[index1];
        
        //double E = pairwiseMinEnergyMatrix.getIntraAndShellE(index1);
        double E = pem.getIntraAndShellE(level, aa1, rc1);

        //for(int level2=level+1; level2<numTreeLevels; level2++){
        for(int level2=0; level2<level; level2++){
            double bestE = Double.POSITIVE_INFINITY;//best pairwise energy

            for( int rot : levelOptions.get(level2).get(conf[level2]) ){
                int index2 = pairwiseMinEnergyMatrix.resOffsets[level2] + rot;
                //bestE = Math.min( bestE, pairwiseMinEnergyMatrix.getPairwiseE(index1,index2) );
                int aa2 = pairwiseMinEnergyMatrix.indicesEMatrixAA[index2];
                int rc2 = pairwiseMinEnergyMatrix.indicesEMatrixRot[index2];
                bestE = Math.min( bestE, pem.getPairwiseE(level,aa1,rc1,level2,aa2,rc2) );
            }      

            E += bestE;
        }

        return E;
    }
    
    
    double quickFScore(int conf[]){
        //get a quick min-sum-min-based lower bound for a node's conformational energies
        double f = 0;
            
        for(int level=0; level<numTreeLevels; level++){

            double bestE = Double.POSITIVE_INFINITY;//best energy (w.r.t. rotamer choices)
            //for the current level, including interactions with all
            //subsequent residues

            for( int rot : levelOptions.get(level).get(conf[level]) ){
                int index1 = pairwiseMinEnergyMatrix.resOffsets[level] + rot;
                bestE = Math.min(bestE, rotE(level,index1,conf) );
            }                

            f += bestE;
        }
        
        //TEMPORARILY NOT DOING LAZY F SCORE
        //THIS MEANS ALWAYS USING "CORRECT" SCORE
        //(CONVERGES TO EXACT ENERGY AS WE DECIDE ALL RES)
        //FOR BOTH LB AND UB
        if(cetm!=null){
            double contContrib = contGScore(conf);
            f += contContrib;
            
            
            //DEBUG!!!!
            /*if (f<-5000){
                
                System.out.println("ERROR: f="+f+" contContrib="+contContrib+" Conf: ");
                for(int c : conf)
                    System.out.print(c + " ");
                System.out.println();
                KSParser.outputObject(this, "CONFTREE5000.dat");
                System.exit(1);
                
                int qqq = 0;
                if(getFirstSplittableLevel(conf)==numTreeLevels){//fully expanded and still so negative...hmm
                    int aaa = 0;
                }
            }*/
            //DEBUG!!!
            
        }

        //DEBUG!!!
        /*
        for(int a : conf)
            System.out.print(a+" ");
        
        //compare to MPLP?  this works in special case with full assignments at ress
        int[] confm = new int[numTreeLevels];
        for(int level=0; level<numTreeLevels; level++){
            if(conf[level]<numSimpleOptions[level])
                confm[level] = conf[level];
            else
                confm[level] = -1;
        }
        
        MSMPLP MPLPheuristic = new MSMPLP(numTreeLevels,numSimpleOptions,pairwiseMinEnergyMatrix);
        double mplpe = MPLPheuristic.optimizeEMPLP(confm, EnvironmentVars.MPLP_iterations);
        double pblpe = fScorePBLP(conf,true);//TEMPORARILY INCLUDING CONT FOR "QUICK" PBLP
        
        
        //also want cont component of g-score (though not part of quick score...)
        double contg = contGScore(conf);
        
        System.out.println(" Quick f-score: "+f+" LP-based quick f-score: "+pblpe+" MPLP: "+mplpe
                +" cont component of g-score: "+contg);
        */
        //DEBUG!!!
        
        return f;
    }
    
    
    double rotE(int level, int index1, int[] conf){
        //contribution to quickFScore for given rotamer (reduced, long index) at given level
        double E = pairwiseMinEnergyMatrix.getIntraAndShellE(index1);

        //for(int level2=level+1; level2<numTreeLevels; level2++){
        for(int level2=0; level2<level; level2++){
            double bestE = Double.POSITIVE_INFINITY;//best pairwise energy

            for( int rot : levelOptions.get(level2).get(conf[level2]) ){
                int index2 = pairwiseMinEnergyMatrix.resOffsets[level2] + rot;
                bestE = Math.min( bestE, pairwiseMinEnergyMatrix.getPairwiseE(index1,index2) );
            }      

            E += bestE;
        }

        return E;
    }
    
    
    double contGScore(int[] conf){
        
        numFS++;//this is basically counting continuous g-score evaluations...
        
        int AAind[][] = new int[numTreeLevels][], rot[][] = new int[numTreeLevels][];
        expandConfToAARot(conf,AAind,rot);
        
        //assuming for now that we expand in usual order!!!
        int splitRes[] = new int[getFirstSplittableLevel(conf)];
        for(int s=0; s<splitRes.length; s++)
            splitRes[s] = s;
        
        PartialMinConstrainer pmc = new PartialMinConstrainer(pemUnreduced,null,splitFlags,tripleFlags,AAind,rot);
        pmc.contOnly = true;
        pmc.cetm = cetm;        
        
        double score = pmc.minTupleE(splitRes);
        
        //The score should always be positive; negative scores indicate numerical error
        //this is unlikely to be significant except in cases with very high variation across voxel,
        //indicating a steric clash
        //so in these cases having a 0 contGScore, and letting the large scalar component keep us out of
        //the voxel, will suffice
        if(score<0)
            score = 0;
        
        return score;
    }
    
    
    double fScorePBLP(int conf[], boolean quick){
        //lower bound for node's conformational energies built from bounds on subsets of the energy
        //(like total pair energies) by linear programming
        int AAind[][] = new int[numTreeLevels][], rot[][] = new int[numTreeLevels][];
        expandConfToAARot(conf,AAind,rot);
        
        PartialMinConstrainer pmc;
        if(quick){
            pmc = new PartialMinConstrainer(pemUnreduced,null,splitFlags,tripleFlags,AAind,rot);

            pmc.cetm = cetm;//DEBUG!!!
        }
        else
            throw new RuntimeException("NOT SUPPORTED YET");
            //pmc = new PartialMinConstrainer(pemUnreduced,cetm);
        
        double ans = pmc.getBound();
        return ans;
    }
    
    
    void expandConfToAARot(int conf[], int AAind[][], int rot[][]){
        //Given conf, in the assignment list form used in tree nodes,
        //expand to lists of AA indices and rotamers
        //so residue res has assignment conf[res], which is expanded so that for each j
        //(AAind[res][j],rot[res][j]) is a rotamer in that assignment
        for(int res=0; res<numTreeLevels; res++){
            int numRot = levelOptions.get(res).get(conf[res]).size();
            AAind[res] = new int[numRot];
            rot[res] = new int[numRot];
            for(int j=0; j<numRot; j++){
                int index = pairwiseMinEnergyMatrix.resOffsets[res] + levelOptions.get(res).get(conf[res]).get(j);
                AAind[res][j] = pairwiseMinEnergyMatrix.indicesEMatrixAA[index];
                rot[res][j] = pairwiseMinEnergyMatrix.indicesEMatrixRot[index];
            }
        }
    }
    
    
    
    ArrayList<ConfTreeNode> getChildren(ConfTreeNode expNode, ConfSpaceSplit split){
        //make children of expNode using specified split
        //and assigning lower bounds with quick f-scores
        ArrayList<ConfTreeNode> ans = new ArrayList<>();
        
        for(int childOption : split.splitOptions){
            //we consider the child with childOption at res
            //we assume split is specially made for expNode, so all simple options in childOption
            //are available at expNode
            int childConf[] = expNode.conf.clone();
            childConf[split.res] = childOption;
            
            if(!canPrune(childConf)){
                ans.add(makeNode(childConf));
            }
        }
        
        return ans;
    }
    
    
    boolean canPrune(int conf[]){
        //if we have split flags, we can try to prune new nodes
        //basically we'll prune if there are no unpruned pairs at some pair of positions in conf
        if(splitFlags!=null){
            for(int res=0; res<numTreeLevels; res++){
                for(int res2=0; res2<res; res2++){
                    boolean hasConf = false;
                    
                    for(int rot : levelOptions.get(res).get(conf[res])){
                        
                        int AA1 = pairwiseMinEnergyMatrix.indicesEMatrixAA[pairwiseMinEnergyMatrix.resOffsets[res]+rot];
                        int rc1 = pairwiseMinEnergyMatrix.indicesEMatrixRot[pairwiseMinEnergyMatrix.resOffsets[res]+rot];
                        
                        for(int rot2 : levelOptions.get(res2).get(conf[res2]) ){
                            
                            int AA2 = pairwiseMinEnergyMatrix.indicesEMatrixAA[pairwiseMinEnergyMatrix.resOffsets[res2]+rot2];
                            int rc2 = pairwiseMinEnergyMatrix.indicesEMatrixRot[pairwiseMinEnergyMatrix.resOffsets[res2]+rot2];

                            if(!splitFlags[res][AA1][rc1][res2][AA2][rc2]){
                                hasConf = true;
                                break;
                            }
                        }
                        
                        if(hasConf)
                            break;
                    }
                    
                    if(!hasConf)//no conf possible for this res pair!
                        return true;
                }
            }
        }
        
        return false;//can't prune
    }
    
    
    ConfTreeNode makeNode(int conf[]){
        //make node with quick F-score
        ConfTreeNode ans = new ConfTreeNode();
        ans.conf = conf;
        ans.LB = quickFScore(conf);
        return ans;
    }
    
    
    
    boolean updateUB(ConfTreeNode expNode, int[] startingConf){
        //Get an upper-bound on the node by a little FASTER run, generating UBConf
        //store UBConf and UB in expNode
        //we'll start with the starting conf (likely from a parent) if provided
        //return true unless no valid conf is possible...then false
        
        int[] UBConf = startingConf.clone();
        //ok first get rid of anything in startingConf not in expNode's conf space,
        //replace with the lowest-intra+shell-E conf
        for(int level=0; level<numTreeLevels; level++){
            if( ! levelOptions.get(level).get(expNode.conf[level]).contains( startingConf[level] ) ){
                
                double bestISE = Double.POSITIVE_INFINITY;
                int bestRC = levelOptions.get(level).get(expNode.conf[level]).get(0);
                for(int rc : levelOptions.get(level).get(expNode.conf[level]) ){
                    double ise = pairwiseMinEnergyMatrix.getIntraAndShellE(level, rc);
                    if( ise < bestISE){
                        bestISE = ise;
                        bestRC = rc;
                    }
                }
                
                UBConf[level] = bestRC;
            }
        }
        
        if(splitFlags!=null){
            if(!makeSureUnpruned(UBConf))
                return false;
        }
        
        double curE = quickFScore(UBConf);
        
        
        if(rigidMtx!=null){//special rigid matrix provided
            //let's refine with rigid, see if we can improve on UBConf
            int UBConf2[] = UBConf.clone();
            double curERigid = fScoreRigid(UBConf,rigidMtx);
            
            boolean done = false;

            int numRounds = 0;//DEBUG!!

            while(!done){

                done = true;

                for(int level=0; level<numTreeLevels; level++){

                    int testConf[] = UBConf2.clone();

                    for(int rc : levelOptions.get(level).get(expNode.conf[level]) ){
                        testConf[level] = rc;
                        
                        if(!canPrune(testConf)){//pruned conf unlikely to be good UB
                        
                            double testE = fScoreRigid(testConf,rigidMtx);
                            if( testE < curERigid){
                                curERigid = testE;
                                UBConf2[level] = rc;
                                done = false;
                            }
                        }
                    }
                }

                numRounds++;
                //if(numRounds==10)
                //    System.out.println("Got to 10 rounds in updateUB...");
            }
            
            //see if low-rigid-E solution is better in continuous space
            double curE2 = quickFScore(UBConf2);
            
            if(curE2<curE){
                UBConf = UBConf2;
                curE = curE2;
            }
        }
        

        boolean refineUB = false;
        
        if(refineUB){
            //now let's do one-res conf changes until convergence.  
            //we use the fact that quickFScore(fully assigned conf) = energy, FOR CURRENT RIGID SCHEME
            //(CAN SPEED THIS UP BY DEDICATED FULLY-ASSIGNED CONF E EVALUATOR IF NEEDED)
            boolean done = false;

            int numRounds = 0;//DEBUG!!

            while(!done){

                done = true;

                for(int level=0; level<numTreeLevels; level++){

                    int testConf[] = UBConf.clone();

                    for(int rc : levelOptions.get(level).get(expNode.conf[level]) ){
                        testConf[level] = rc;
                        
                        if(!canPrune(testConf)){//pruned conf unlikely to be good UB
                        
                            double testE = quickFScore(testConf);
                            if( testE < curE){
                                curE = testE;
                                UBConf[level] = rc;
                                done = false;
                            }
                        }
                    }
                }

                numRounds++;
                if(numRounds==10)
                    System.out.println("Got to 10 rounds in updateUB...");
            }
        }
        
        expNode.UBConf = UBConf;
        expNode.UB = curE;
        
        return true;
    }
    
    
    boolean makeSureUnpruned(int conf[]){
        //make changes to (fully-assigned) conf, if needed, until it has no pruned pairs
        //if impossible return false (indicates all confs pruned)
        
        return makeSureUnprunedHelper(conf,0);
    }
    
    boolean makeSureUnprunedHelper(int[] conf, int level){
        //assuming there are no pruned pairs in conf among residues < level
        //see if we can fill in the rest of the conformation without hitting any pruned pairs
        //at each level, we only change RCs in conf if necessary to avoid pruning
        
        if(level==numTreeLevels)//made it here with no pruning!
            return true;
        
        int origOpt = conf[level];

        if(confOKUpToLevel(conf,level)){
            if(makeSureUnprunedHelper(conf,level+1))
                return true;
        }
            
        for(int opt=0; opt<numSimpleOptions[level]; opt++){
            if(opt!=origOpt){
                conf[level] = opt;

                if(confOKUpToLevel(conf,level)){
                    if(makeSureUnprunedHelper(conf,level+1))
                        return true;    
                }
            }
        }
        
        conf[level] = origOpt;//found nothing...return conf to original state
        //(if level>0, we'll now backtrack and try making changes at lower levels)
        return false;
    }
    
    
    boolean confOKUpToLevel(int[] conf, int level){
        //Does (fully-defined) conf have any pruned pairs among residues up to specified level?
        
        for(int level2=0; level2<level; level2++){

            int AA1 = pairwiseMinEnergyMatrix.indicesEMatrixAA[pairwiseMinEnergyMatrix.resOffsets[level]+conf[level]];
            int rc1 = pairwiseMinEnergyMatrix.indicesEMatrixRot[pairwiseMinEnergyMatrix.resOffsets[level]+conf[level]];
            int AA2 = pairwiseMinEnergyMatrix.indicesEMatrixAA[pairwiseMinEnergyMatrix.resOffsets[level2]+conf[level2]];
            int rc2 = pairwiseMinEnergyMatrix.indicesEMatrixRot[pairwiseMinEnergyMatrix.resOffsets[level2]+conf[level2]];

            if(splitFlags[level][AA1][rc1][level2][AA2][rc2])
                return false;
        }
        
        return true;
    }
    
    

}
