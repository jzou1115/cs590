
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;
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
//	SeqTree.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
//This is a tree of sequences.  Meant for optimizations in sequences space
//where both the objective function and any constraints are functions of sequence
//defined as linear functions of the GMEC energies for different states. 
//States described in the constituent ConfTrees (e.g., different ligands bound).  
//for example, minimize binding enthalpy to ligand 1 while constraining binding enthalpy to ligand 2
public class SeqTree implements Serializable {
    
    //information on optimization problem (variables, constraints, objective function)
    
    int numTreeLevels;//number of residues with sequence changes
    
    GMECLinFunc objFcn;//we are minimizing objFcn...
    GMECLinFunc[] constraints;//with the constraints that constr all <= 0
    
    int numAATypes[];//number of allowed AA types at each mutable residue
    ArrayList<ArrayList<ArrayList<Integer>>> AATypeOptions;
    //SeqTreeNode.conf assigns each level an index in AATypeOptions.get(level)
    //At this index is a list of available AA types (rotamer library numbering)
    //allowing conf to specify the sequence space.  Similar to ConfTree.levelOptions
    //NOTE: CONF IS NOT A LIST OF AA TYPE INDICES (AS IN ROTAMER LIBRARY...USE AATYPEOPTIONS TO CONVERT TO THAT)!!!
    
    
    //information on states
    
    int numStates;//how many states there are
    //they have to have the same mutable residues & AA options,
    //though the residues involved may be otherwise different
    
    int stateNumRes[];//number of flexible residues in each state
    ArrayList<ArrayList<Integer>> mutable2StateResNums;
    //mutable2StateResNum.get(state) maps levels in this tree to flexible residues for state
    //(not necessarily an onto mapping)
    
    ArrayList<ArrayList<ArrayList<int[]>>> stateRCs;//stateRCs.get(state).get(res) for a given state
    //lists (AA index, RC index) pairs avaible for this res (numbered among flexible res for state)
    ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> stateRCsAtAA;
    //stateRCsAtAA.get(state).get(res).get(AAindex)
    //lists available RC indices for that AA type.  Derived from stateRCs
    
    PairwiseEnergyMatrix stateEMatrix[];//energy matrix for each state
    

    CETMatrix stateCETM[];//If not null, we use these to model cont flex
    //not currently used for negative-coefficient states for partially-defined seq though 
    //(which is valid because discr GMEC>cont GMEC for each seq)

    
    //now, for negative-coefficient states we'll currently need to use discrete rots to bound
    //for partially defined sequences
    //we'll have different pruning & E-matrix, and maybe a different rot lib
    ArrayList<ArrayList<ArrayList<int[]>>> negStateRCs;
    ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> negStateRCsAtAA;
    PairwiseEnergyMatrix negStateEMatrix[];
    
    
    int numMaxMut;//number of mutations allowed away from wtSeq (-1 means no cap)
    int wtSeq[];
    
    
    RotamerLibrary rl;//a rotamer library for the mutations
    //will be taken from first state
    
    
    //The queue and statistics on it
    
    PriorityBlockingQueue<SeqTreeNode> seqExpansion = new PriorityBlockingQueue<>();
    //this doesn't need to be a queue because we advance based on the stateTrees' queues

    int numExpanded = 0;
    int numFS = 0;
    int numPruned = 0;//nodes are pruned if we know they contain no confs w/i the constraints
    int numSeqDefExpansions = 0;//number of expanded nodes creating fully defined sequences
    int numSeqDefNodes = 0;//how many nodes with fully defined sequence are in the queue
    int numSeqsReturned = 0;
    int stateGMECsForPruning = 0;//how many state GMECs have been calculated for nodes that are pruned
    
    //ordering controls
    static boolean considerConfBest = true;//consider only conf splits that are optimal for their tree

    static boolean fixedOrder = true;
    //this is the simplest ordering scheme:
    //first expand the sequence, one res at a time, getting lower bounds using calcLSSubstQuick
    //(which uses (type-dependent) pruning of RCs but no splitting of confs w/i sequence)
    //then we start expanding confs at each state, favoring states w/ more res still to expand
    //during this phase (fully-defined seq), we can use ConfTree's lowest nodes' LBs and UBs to bound
    //state GMEC energies
    
    
    static boolean lazyUB = true;//only call updateUB on lowest node in tree when we need to do a seqNode
    //LB, instead of calling for each new node
    
    //To do SeqTree A*:
    //we're minimizing function of seqExpansion, which is lower-bounded in seqExpansion
    //at returning point, f score must be actual obj function and must be w/i constraints
    //and this is for lowest-scoring node
    //so we should always be expanding lowest-scoring seq node
    //we can do this by expanding the lowest conf node in any of its states
    //we have a judgment call: consider only best split (in terms of converging to GMEC) for each state conf tree
    //at our current seq node, and pick the one that advances seqExpansion best,
    //OR looks at all split options and pick the one that advances seqExpansion best 
    //(for this could even consider ordering tree by UB)
    //MIGHT TRY BOTH BUT START BY LOOKING AT BESTSPLITS ONLY
    
    
    public SeqTree(RotamerSearch[] stateRS, RotamerSearch negStateRS[], int stateStrandMut[][][], ParamSet sParams) {
        //Given RotamerSearch and strandMut objects set up for each state (with energy matrix loaded,
        //as if to do a GMEC search for that state)
        //set up a sequence tree and initialize the sequence expansion with a totally unassigned node
        //We also read general search (rather than state-specific) info directly from sParams
        //negStateRS are for use in bounding functions at partially defined sequence nodes, for states with coeffs<0
        //they can point to stateRS if appropriate
        
        //INPUT INTERFACE: PROVIDE SYSTEM.CFG AND DEE.CFG FOR EACH STATE JUST AS IF DOING GMEC SEARCH
        //THEN PROVIDE Multistate.cfg WITH GENERAL INFO, WHICH WILL BE READ HERE
        
        numTreeLevels = (new Integer((String)sParams.getValue("NUMMUTRES"))).intValue();
        int numConstr = (new Integer((String)sParams.getValue("NUMCONSTR"))).intValue();
        numStates = (new Integer((String)sParams.getValue("NUMSTATES"))).intValue();
        
        objFcn = new GMECLinFunc((String)sParams.getValue("OBJFCN"));
        
        constraints = new GMECLinFunc[numConstr];
        for(int constr=0; constr<numConstr; constr++)
            constraints[constr] = new GMECLinFunc((String)sParams.getValue("CONSTR"+constr));
        
        
        mutable2StateResNums = new ArrayList<>();
        //these will be listed directly in Multistate.cfg under "STATEMUTRES0" etc.
        for(int state=0; state<numStates; state++){
            ArrayList<Integer> m2s = new ArrayList<>();
            String stateMutRes = (String)sParams.getValue("STATEMUTRES"+state);
            StringTokenizer st = new StringTokenizer(stateMutRes);
            
            while(st.hasMoreTokens())
                m2s.add(Integer.valueOf(st.nextToken()));
            
            if(m2s.size()!=numTreeLevels){
                throw new RuntimeException("ERROR: SeqTree has "+numTreeLevels+" mutable positions "
                        +" but "+m2s.size()+" are listed for state "+state);
            }
            
            mutable2StateResNums.add(m2s);
        }
        
        //we can collect info on allowed sequences from any state (we choose state 0)
        //but we'll check that they all agree
        parseAllowedSeq(stateRS[0],stateStrandMut[0],0,false);
           
        //fill in information on states
        stateNumRes = new int[numStates];
        stateEMatrix = new PairwiseEnergyMatrix[numStates];
        stateRCs = new ArrayList<>();
        stateRCsAtAA = new ArrayList<>();
        
        negStateEMatrix = new PairwiseEnergyMatrix[numStates];
        negStateRCs = new ArrayList<>();
        negStateRCsAtAA = new ArrayList<>();
        
        stateCETM = new CETMatrix[numStates];
        
        for(int state=0; state<numStates; state++){
            //check consistency of sequence options
            parseAllowedSeq(stateRS[state],stateStrandMut[state],state,true);
            
            stateNumRes[state] = stateRS[state].numberMutable;
            stateEMatrix[state] = stateRS[state].getMinMatrix();
            
            stateCETM[state] = stateRS[state].cetm;
            
            makeRCLists(stateRS[state],stateNumRes[state],stateRCs,stateRCsAtAA);
            
            if(stateRS[state]==negStateRS[state]){
                //no distinct pruning, etc. for negative
                //just link things
                negStateRCs.add(stateRCs.get(state));
                negStateRCsAtAA.add(stateRCsAtAA.get(state));
                negStateEMatrix[state] = stateEMatrix[state];
            }
            else {
                //load distinct negative info
                negStateEMatrix[state] = negStateRS[state].getMinMatrix();
                makeRCLists(negStateRS[state],stateNumRes[state],negStateRCs,negStateRCsAtAA);
            }     
            
        }
        
        
        rl = EnvironmentVars.aaRotLib;//ASSUMING FOR NOW AMINO ACID MUTATIONS

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)
        numMaxMut = (new Integer((String)sParams.getValue("NUMMAXMUT","-1"))).intValue();
        if(numMaxMut>-1){
            wtSeq = parseWTSeq((String)sParams.getValue("WTSEQ"));
        }
        
        
        //set up empty expansion queue, like in ConfTree A* for run1
        
        //start by just putting in all undecided conformation
        int conf[] = new int[numTreeLevels];
        for(int level=0; level<numTreeLevels; level++)
            conf[level] = numAATypes[level];

        SeqTreeNode seqNode = new SeqTreeNode();
        seqNode.conf = conf;

        seqNode.statePrunedRot = new DEEGoldstein[numStates];
        seqNode.statePrunedPairs = new DEEGoldsteinPairs[numStates];
        seqNode.stateTrees = new ConfTree[numStates];//let's hold off on these for now...
        //build when sequence fully defined
        
        seqNode.negStatePrunedRot = new DEEGoldstein[numStates];
        seqNode.negStatePrunedPairs = new DEEGoldsteinPairs[numStates];

        //build pruning structures...
        for(int state=0; state<numStates; state++){
            
            seqNode.setupPruning(stateRS[state],false,state,stateStrandMut);
            
            if(stateRS[state]==negStateRS[state]){
                seqNode.negStatePrunedRot[state] = seqNode.statePrunedRot[state];
                seqNode.negStatePrunedPairs[state] = seqNode.statePrunedPairs[state];
            }
            else
                seqNode.setupPruning(negStateRS[state],true,state,stateStrandMut);
        }

        seqNode.LB = calcLB(seqNode,objFcn);
        seqExpansion.put(seqNode);
    }
    
    
    int[] parseWTSeq(String seq){
        //we'll convert this sequence, specified like ARG LYS...,
        //to AA type indices based on the first state's rot lib
        StringTokenizer st = new StringTokenizer(seq);
        if(st.countTokens()!=numTreeLevels)
            throw new RuntimeException("ERROR: wrong number of residues in WT seq: "+seq);
        
        int wt[] = new int[numTreeLevels];
        
        for(int level=0; level<numTreeLevels; level++)
            wt[level] = rl.getAARotamerIndex(st.nextToken());
        
        return wt;
    }
    
    
    void makeRCLists( RotamerSearch rs, int numRes, ArrayList<ArrayList<ArrayList<int[]>>>
            RCList, ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> RCAtAAList ){
        //get lists of RCs from the given RotamerSearch, with the given number of residues
        //add to the RC lists (intended to be either stateRCs, stateRCsAtAA or the neg versions)
        
        //preallocate stateRCs and stateRCsAtAA to have a list for every (state, res)
        ArrayList<ArrayList<int[]>> stateRCList = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> stateRCListAA = new ArrayList<>();
        
        for(int res=0; res<numRes; res++){
            stateRCList.add(new ArrayList<int[]>());
            stateRCListAA.add(new ArrayList<ArrayList<Integer>>());
        }

        //now let's record all the RCs in the rotamer search
        //eliminatedRotAtRes is a convenient place to find them
        for(RotInfo<Boolean> ri : rs.eliminatedRotAtRes){
            stateRCList.get(ri.curPos).add(new int[] {ri.curAA,ri.curRot});

            while(stateRCListAA.get(ri.curPos).size()<=ri.curAA)//allocate list of RCs for AA type if needed
                stateRCListAA.get(ri.curPos).add(new ArrayList<Integer>());

            stateRCListAA.get(ri.curPos).get(ri.curAA).add(ri.curRot);
        }
        
        RCList.add(stateRCList);
        RCAtAAList.add(stateRCListAA);
    }
    
    
    void parseAllowedSeq(RotamerSearch rs, int[][] strandMut, int state, boolean check){
        //figure out numAATypes and AATypeOptions, using the information for one of the states
        //if check then instead of filling in the information, we check to make sure
        //it matches what's already there (used to ensure all states agree; we assume state 0 used to fill in first)
        
        if(!check){
            numAATypes = new int[numTreeLevels];
            AATypeOptions = new ArrayList<>();
        }
        
        for(int level=0; level<numTreeLevels; level++){
            
            ArrayList<ArrayList<Integer>> levelAAOpt = new ArrayList<>();
            
            int res = mutable2StateResNums.get(state).get(level);//residue number for state 0 for this level
            int strandNum = rs.mutRes2Strand[res];
            int strandResNum = strandMut[strandNum][rs.mutRes2StrandMutIndex[res]];
            
            int AATypeCount = rs.strandRot[strandNum].getNumAllowable(strandResNum);
            if(check){
                if(AATypeCount!=numAATypes[level]){
                    throw new RuntimeException("ERROR: Level "+level+" has "+numAATypes[level]
                            +" sequence options for first state but "+AATypeCount+" for state "+state);
                }
            }
            else
                numAATypes[level] = AATypeCount;
            
            ArrayList<Integer> allAllowed = new ArrayList<>();
            
            for(int aa=0; aa<numAATypes[level]; aa++){
                ArrayList<Integer> opt = new ArrayList<>();//this is the assignment option for the current AA type
                int curAA = rs.strandRot[strandNum].getIndexOfNthAllowable(strandResNum,aa);
                opt.add(curAA);
                allAllowed.add(curAA);
                levelAAOpt.add( opt );
                
                if(check){
                    if(curAA!=AATypeOptions.get(level).get(aa).get(0)){
                        throw new RuntimeException("ERROR: Level "+level+" AA type "+aa+" is "+
                                AATypeOptions.get(level).get(aa).get(0)+" for first state but "
                                +curAA+" for state "+state);
                    }
                }
            }
            
            levelAAOpt.add(allAllowed);//add an unassigned-AA-type option
            
            if(!check)
                AATypeOptions.add(levelAAOpt);
        }

    }
    
    
    //compute the best sequence remaining in the tree (i.e., in the union of the nodes's seq spaces)
    public int[] doAStar ( /*int numMaxChanges, int nodesDefault[], 
			String strandDefault[][], */){

        /*
         * While seq not fully defined
         * 	For the current minimum node in the queue
         * 		Consider all splits at any of the trees; pick best and apply it
         *              Recalculate score at corresponding SeqTreeNode
         *              If the best imposes an AA type at a previously undefined res
         *                      Split SeqTreeNode accordingly, calculating scores for child nodes
         * 			Propagate split to other trees (CAN BE LAZY)
         */
        while (true) {	

            SeqTreeNode seqNode = seqExpansion.poll();
            
            
            if (seqNode==null){//the queue is empty
                int emptyConf[] = new int[numTreeLevels];
                Arrays.fill(emptyConf,-1);
                System.out.println("No more unpruned nodes in SeqTree...terminating A*");
                return emptyConf;//so return a sequence of -1's to flag the stop of the search
            }

            if( seqNode.statePrunedPairs==null && (!isSeqFullyDefined(seqNode)) ){
                //need to decompress pruned pairs so we can use them
                //for now this is only needed if not fully defined
                seqNode.statePrunedPairs = new DEEGoldsteinPairs[numStates];
                seqNode.negStatePrunedPairs = new DEEGoldsteinPairs[numStates];
                
                for(int state=0; state<numStates; state++){
                    
                    DEEGoldsteinPairs uncPrunedPairs = uncompressPrunedPairs(seqNode.npp[state],
                            seqNode.statePrunedRot[state].eliminatedRotAtPos,
                            state, null);
                    
                    
                    /*
                    for(int res=0; res<stateNumRes[state]; res++){
                        for(int res2=0; res2<res; res2++){
                            for(int[] rc : stateRCs.get(state).get(res)){
                                for(int[] rc2 : stateRCs.get(state).get(res2)){
                                    if(uncPrunedPairs.splitFlags[res][rc[0]][rc[1]][res2][rc2[0]][rc2[1]]
                                            != seqNode.statePrunedPairs[state].splitFlags[res][rc[0]][rc[1]][res2][rc2[0]][rc2[1]]){
                                        System.out.println("problem...");
                                    }
                                }
                            }
                        }
                    }*/
                    
                    
                    seqNode.statePrunedPairs[state] = uncPrunedPairs;
                    seqNode.statePrunedRot[state].splitFlags = uncPrunedPairs.splitFlags;
                    
                    if(seqNode.npp[state]==seqNode.negNPP[state]){
                        seqNode.negStatePrunedPairs[state] = uncPrunedPairs;
                        //and negStatePrunedRot will be the same as statePrunedRot
                    }
                    else {
                        seqNode.negStatePrunedPairs[state] = uncompressPrunedPairs(seqNode.negNPP[state],
                            seqNode.negStatePrunedRot[state].eliminatedRotAtPos,
                            state, null);
                        
                        seqNode.negStatePrunedRot[state].splitFlags = seqNode.negStatePrunedPairs[state].splitFlags;
                    }
                }
            }
            

            //ok we now need to decide on a split in one of the conf states
            //and of course this goes with one of the trees
            ConfSpaceSplit split = pickBestSplit(seqNode);

            if(split==null){//no splits available for any conf trees
                //thus, the sequence is fully defined, and each ConfTree's lowest node
                //is the GMEC for this sequence, so seqNode's score is the actual 
                //function value for the sequence (rather than a non-tight lower bound)
                //since this is the lowest-scored node, we have found the optimal sequence
                //and should return it
                
                numSeqDefNodes--;
                numSeqsReturned++;
                printBestSeqInfo(seqNode);   
                
                
                
                //DEBUG!!!!!
                //KSParser.outputObject(this, "SEQTREEDONE.dat");
                //KSParser.outputObject(seqNode, "BESTSEQNODE.dat");
                
                
                return seqNode.conf;
            }
            //split can either be a sequence split or not
            else if (split.parentConfTree==null){//indicates a sequence split
                
                boolean fullyDefining = false;
                
                for(ConfTreeNode child : getChildren(seqNode,split)){//child is actually a SeqTreeNode...
                    SeqTreeNode schild = (SeqTreeNode)child;
                    
                    if(!canPrune(schild)){
                        
                        compressPrunedPairs(schild,seqNode);                        
                        
                        seqExpansion.add(schild);
                        if(isSeqFullyDefined(schild)){
                            numSeqDefNodes++;
                            fullyDefining = true;
                        }
                    }
                    else{
                        numPruned++;
                        stateGMECsForPruning += countStateGMECs(schild);
                    }
                }
                
                numExpanded++;
                if(fullyDefining)
                    numSeqDefExpansions++;
                    
                
                if( numExpanded%200/*1000*/ == 0 ){
                    System.out.println(numExpanded+" expanded, of which "+numSeqDefExpansions
                            +" fully defined sequences; "+seqExpansion.size()+" nodes, of which "
                            +numSeqDefNodes+" are fully defined; "+
                            numPruned+" pruned, lowest bound="+seqExpansion.peek().LB);
                    
                    
                    
                    //DEBUG!!  Just to pause
                    /*if(numExpanded==15000){
                        try{System.in.read();}
                        catch(Exception e){}
                    }*/
                    
                    
                }
            }
            else {
                ConfTree parentTree = split.parentConfTree;
                ConfTreeNode nodeToSplit = parentTree.curExpansion.poll();
                
                ArrayList<ConfTreeNode> children = parentTree.getChildren(nodeToSplit,split);
                            
                for(ConfTreeNode child : children){

                    if(child!=null){
                        boolean addNode = true;

                        //NOT YET DOING MAX NUMBER OF MUTATIONS
                        /*
                        if(!singleSeq){
                            if(countMutants(child,nodesDefault)>numMaxChanges)
                                addNode = false;
                        }*/

                        if( addNode ){
                            
                            if(!lazyUB){
                                //we'll be needing upper as well as lower bounds
                                if(!parentTree.updateUB(child,nodeToSplit.UBConf))
                                    continue;//node has inevitable clash...don't add it to the expansion

                                int state = getStateNum(seqNode,parentTree);
                                seqNode.stateUB[state] = Math.min(seqNode.stateUB[state],child.UB);
                                //The minimum GMEC upper bound for any node for this sequence and state
                                //is a valid upper bound on the sequence GMEC in this state,
                                //because one of the nodes must contain the GMEC
                            }
                            else{
                                child.UB = Double.NaN;//flag we didn't calculate the UB yet...
                                child.UBConf = nodeToSplit.UBConf;//to use later
                            }
                            
                            parentTree.curExpansion.put(child);
                        }
                    }
                }
                                
                seqNode.LB = calcLB(seqNode,objFcn);
                if(!canPrune(seqNode)){
                    
                    if(seqNode.statePrunedPairs!=null)//we have uncompressed pruning info
                        //so we need to update the compressed version
                        compressPrunedPairs(seqNode,null);
                    
                    seqExpansion.put(seqNode);//put back in with new score
                }
                else{
                    numSeqDefNodes--;
                    numPruned++;
                    stateGMECsForPruning += countStateGMECs(seqNode);
                }
            }
        }
    }
    
    
    void printBestSeqInfo(SeqTreeNode seqNode){
        //About to return the given fully assigned sequence from A*
        //provide information
        System.out.println("SeqTree: A* returning conformation; lower bound = "+seqNode.LB+" nodes expanded: "+numExpanded+" FS terms evaluated: "+numFS);
        
        System.out.print("Sequence: ");

        for(int level=0; level<numTreeLevels; level++)
            System.out.print(AATypeOptions.get(level).get(seqNode.conf[level]).get(0)+" ");
        System.out.println();
        
        //let's output in AA name form...
        for(int level=0; level<numTreeLevels; level++)
            System.out.print( rl.getAAName(AATypeOptions.get(level).get(seqNode.conf[level]).get(0))+" ");
        System.out.println();
        
        //provide state GMECs, specified as rotamers (AA types all the same of course)
        for(int state=0; state<numStates; state++){
            System.out.print("State "+state);
            
            if(seqNode.stateTrees[state]==null) {
                System.out.println(" has an unavoidable clash.");
            }
            else {
                System.out.print(" RCs: ");
                int conf[] = seqNode.stateTrees[state].curExpansion.peek().conf;
                ReducedEnergyMatrix rem = seqNode.stateTrees[state].pairwiseMinEnergyMatrix;
                for(int res=0; res<stateNumRes[state]; res++)
                    System.out.print( rem.indicesEMatrixRot[rem.resOffsets[res]+conf[res]]+" ");

                System.out.println( "Energy: "+
                        (seqNode.stateTrees[state].curExpansion.peek().LB+stateEMatrix[state].getShellShellE()) );
            }
        }
        
        System.out.println();
        System.out.println(numExpanded+" expanded, of which "+numSeqDefExpansions
                            +" fully defined sequences; "+seqExpansion.size()+" nodes in tree, of which "
                            +numSeqDefNodes+" are fully defined; "+
                            numPruned+" pruned.");
        

        int stateGMECsRet = numSeqsReturned*numStates;
        int stateGMECsInTree = countGMECsInTree();
        int totGMECsCalcd = stateGMECsRet + stateGMECsInTree + stateGMECsForPruning;
        System.out.println(totGMECsCalcd+" state GMECs calculated: "+stateGMECsRet+" returned, "+stateGMECsInTree
                +" in tree, "+stateGMECsForPruning+" for pruned sequences.");
    }
    
    
    int countGMECsInTree(){
        //count how many state GMECs have been calculated and are in nodes in the tree--
        //that is, how many ConfTrees at sequence nodes have been expanded to the points
        //that their lowest-bound node is fully expanded (and thus is their GMEC)
        //for comparison, regular A* needs to calculate the GMEC for every (non-clashing)
        //(sequence,state) pair
        //Note this only counts state GMECs currently in the tree
        int count = 0;
        
        for(SeqTreeNode seqNode : seqExpansion){
            count += countStateGMECs(seqNode);
        }
        
        return count;
    }
    
    int countStateGMECs(SeqTreeNode seqNode){
        //how many states for this node have GMECs calculated?
        int count = 0;
        
        for(ConfTree ct : seqNode.stateTrees){
            if(ct!=null){
                int conf[] = ct.curExpansion.peek().conf;
                if(ct.getFirstSplittableLevel(conf)==ct.numTreeLevels)//conf fully defined
                    count++;
            }
        }
        
        return count;
    }
    
    boolean canPrune(SeqTreeNode seqNode){
        //check if the node can be pruned based on the constraints
        //each constraint function must be <=0, so if any constraint function's lower bound
        //over sequences in seqNode is > 0,
        //we can prune seqNode
        
        for(GMECLinFunc constr : constraints){
            if(calcLB(seqNode,constr)>0)
                return true;
        }
        
        if(numMaxMut!=-1){
            //we have a cap on the number of mutations...prune if exceeded
            int mutCount = 0;
            for(int level=0; level<numTreeLevels; level++){
                if(seqNode.conf[level]<numAATypes[level]){//AA type at level is assigned
                    if( AATypeOptions.get(level).get(seqNode.conf[level]).get(0) != wtSeq[level] )//and is different from wtSeq
                        mutCount++;
                }
            }
            
            if(mutCount>numMaxMut)
                return true;
        }
        
        return false;
    }


    ArrayList<ConfTreeNode> getChildren(SeqTreeNode seqNode, ConfSpaceSplit seqSplit){
        //Generate the children of seqNode, split using seqSplit

        ArrayList<ConfTreeNode> ans = new ArrayList<>();
        
        for(int childOption : seqSplit.splitOptions){
            //we consider the child with childOption at res
            //we assume split is specially made for expNode, so all simple options in childOption
            //are available at expNode
            int childConf[] = seqNode.conf.clone();
            childConf[seqSplit.res] = childOption;              
                    
            SeqTreeNode child = new SeqTreeNode();
            child.conf = childConf;
            
            //now split the conf trees and pruning info, and update them
            child.stateTrees = new ConfTree[numStates];
            child.statePrunedRot = new DEEGoldstein[numStates];
            child.statePrunedPairs = new DEEGoldsteinPairs[numStates];
            //child.statePrunedTriples = new DEEGoldsteinTriples[numStates];
            child.stateUB = new double[numStates];
            
            child.negStatePrunedRot = new DEEGoldstein[numStates];
            child.negStatePrunedPairs = new DEEGoldsteinPairs[numStates];
            
            
            for(int state=0; state<numStates; state++){
                //copy over the ConfTree for this state, splitting all the nodes...
                //NOTE: CAN ALSO TRY LAZY SPLITTING OF CONF TREES
                //child.stateTrees[state] = new ConfTree(seqNode.stateTrees[state],seqToConfSplit(seqSplit,seqNode.stateTrees[state]));
                //FOR NOW WILL ONLY MAKE CONFTREE WHEN SEQ FULLY DEFINED; THIS AVOIDS HAVING TO 
                //SPLIT CONF TREES ON SEQ SPLITS
                
                doChildPruning(seqNode,child,state,false,seqSplit,childOption);
                
                if(seqNode.statePrunedRot[state]==seqNode.negStatePrunedRot[state]){
                    child.negStatePrunedRot[state] = child.statePrunedRot[state];
                    child.negStatePrunedPairs[state] = child.statePrunedPairs[state];
                }
                else
                    doChildPruning(seqNode,child,state,true,seqSplit,childOption);
                    //child pruning different (probably rigid) for negative state
                
                if(isSeqFullyDefined(child)){//now that pruning is all done,
                    //switch from pruning mode to ConfTree mode.  Done with seq splits!
                    if(!generateConfTree(child,state)){
                        //the sequence space for this childOption has an unavoidable clash
                        //indicated by previous steric pruning of all RC options at some residue
                        //so the upper and lower bounds for this state and sequence will be considered infinite
                        child.stateTrees[state] = null;
                        child.stateUB[state] = Double.POSITIVE_INFINITY;
                        //clashingSeq = true;
                    }
                }
            }
                        
                        
            //if(!clashingSeq){
            //if(!canPrune(child)){//we handle this in doAStar
                child.LB = calcLB(child,objFcn);
                ans.add( child );
            //}
            //}
        }

        return ans;
    }
    
    
    
    void doChildPruning( SeqTreeNode seqNode, SeqTreeNode child, int state, boolean neg, 
            ConfSpaceSplit seqSplit, int childOption ) {
        //when splitting off child from seqNode, do any additional pruning possible
        //for the specified state.  Do this for negative version if specified
        
        //copy over the pruning info from seqNode
        DEEGoldstein dgp = seqNode.statePrunedRot[state];
        DEEGoldsteinPairs dpp = seqNode.statePrunedPairs[state];
        
        if(neg){
            dgp = seqNode.negStatePrunedRot[state];
            dpp = seqNode.negStatePrunedPairs[state];
        }


        //SLOW COPY...DO W/O STREAMING IF BOTTLENECK
        PrunedRotamers<Boolean> prCopy = null;
        boolean[][][][][][] sfCopy = null;
        try {
            prCopy = (PrunedRotamers<Boolean>)KSParser.deepCopy(dgp.eliminatedRotAtPos);
            sfCopy = (boolean[][][][][][])KSParser.deepCopy(dpp.splitFlags);
        }
        catch(Exception e){
            throw new RuntimeException(e);
        }

        DEEGoldstein dg = new DEEGoldstein( dgp.pairwiseMinEnergyMatrix, dgp.pairwiseMaxEnergyMatrix,
                dgp.numMutable, dgp.strandMut, dgp.curEw-dgp.Ival, dgp.strandRot, prCopy,
                dgp.doMinimize, dgp.indIntMinDEE, dgp.pairIntMinDEE, sfCopy, dgp.useFlags,
                dgp.minimizeBB, dgp.mutRes2Strand, dgp.mutRes2MutIndex, dgp.typeDependent, dgp.doIMinDEE,
                dgp.Ival, dgp.doPerturbations, false );

        //let's not use triples for now (providing triple flags as null)
        boolean tripFlags[][][][][][][][][] = null;

        DEEGoldsteinPairs dp = new DEEGoldsteinPairs( dgp.pairwiseMinEnergyMatrix, dgp.pairwiseMaxEnergyMatrix,
                dgp.numMutable, dgp.strandMut, dgp.curEw-dgp.Ival, dgp.strandRot, prCopy, dpp.resInPair,
                dgp.doMinimize, sfCopy, dpp.useFlags, dpp.magicBullet, dpp.distrDEE, dgp.minimizeBB, 
                dpp.maxScale!=1, dpp.maxScale, dgp.mutRes2Strand, dgp.mutRes2MutIndex, dgp.typeDependent, dgp.doIMinDEE,
                dgp.Ival, tripFlags, dgp.doPerturbations, false );

        //prune the newly excluded AA types from prunedRot
        int splitRes = mutable2StateResNums.get(state).get(seqSplit.res);//split residue numbered for this state

        for( int[] rc : stateRCs.get(state).get(splitRes) ){
            if( ! AATypeOptions.get(seqSplit.res).get(childOption).contains(rc[0]) )//AA type no longer in child's seq space
                dg.eliminatedRotAtPos.set(splitRes, rc[0], rc[1], true);//prune
        }

        //now cycle through rot, pairs, etc until convergence
        do {
            dg.ComputeEliminatedRotConf();
            dp.ComputeEliminatedRotConf();
        } while(dg.numRuns>1 || dp.numRuns>1);//repeat if anything was pruned


        //store these as our new pruning objects
        if(neg) {
            child.negStatePrunedRot[state] = dg;
            child.negStatePrunedPairs[state] = dp;
        }
        else {
            child.statePrunedRot[state] = dg;
            child.statePrunedPairs[state] = dp;
        }
    }
    
    
    boolean generateConfTree(SeqTreeNode seqNode, int state){
        //generate a ConfTree for the given state's conformations within the given SeqTreeNode
        //returning false means the sequence space of seqNode has an unavoidable steric clash
        //as indicated by pruning
        
        ReducedEnergyMatrix eMatrixRed = reduceEnergyMatrix(seqNode,state);
        
        if(eMatrixRed==null)
            return false;
        
        ConfTree tree = new ConfTree(eMatrixRed,stateEMatrix[state]);
        tree.fixedOrder = true;
        if(stateCETM[state]!=null)
            tree.cetm = stateCETM[state];
        
        tree.splitFlags = seqNode.statePrunedPairs[state].splitFlags;
        
        
        //providing rigid matrix in continuous case...
        //DEBUG!!
        if(stateCETM[state]!=null){
            tree.rigidMtx = negStateEMatrix[state];
        }
        
        //initialize expansion queue, with both upper and lower bounds for first node (totally unassigned
        //within the post-pruning conf space for this sequence)
        int conf[] = new int[stateNumRes[state]];
        for(int res=0; res<stateNumRes[state]; res++)
            conf[res] = tree.numSimpleOptions[res];

        ConfTreeNode newNode = tree.makeNode(conf);
        int[] blankConf = new int[tree.numTreeLevels];
        Arrays.fill(blankConf, -1);
        
        if(!tree.updateUB(newNode, blankConf))//all confs clash!
            return false;
            
        
        tree.curExpansion.put(newNode);
        
        seqNode.stateTrees[state] = tree;
        
        seqNode.stateUB[state] = newNode.UB;
        
        return true;
    }
    
    
    ReducedEnergyMatrix reduceEnergyMatrix(SeqTreeNode seqNode, int state){
        //generate a ReducedEnergyMatrix for the given state, based on pruning in seqNode
        //similar to PairwiseEnergyMatrix.reduceMatrix
        //return null if everything pruned at some residue
        
        //we'll include the RCs that haven't been pruned
        PrunedRotamers<Boolean> pr = seqNode.statePrunedRot[state].eliminatedRotAtPos;
        ArrayList<int[]> unprunedRCs = new ArrayList<>();//unpruned RCs as (res, AA index, RC index)
        //everything with the wrong sequence will be pruned already...
        
        
        for(int res=0; res<stateNumRes[state]; res++){
            boolean resOK = false;//set true here when we find an unpruned RC at this res
            
            for(int rc[] : stateRCs.get(state).get(res)){
                if(!pr.get(res,rc[0],rc[1])){
                    unprunedRCs.add(new int[] {res,rc[0],rc[1]});
                    resOK = true;
                }
            }
            
            if(!resOK)
                return null;
        }
        
        
        int totUnpruned = unprunedRCs.size();
        int indicesPos[] = new int[totUnpruned];
        int indicesAA[] = new int[totUnpruned];
        int indicesRot[] = new int[totUnpruned];
        
        for(int q=0; q<totUnpruned; q++){
            int rc[] = unprunedRCs.get(q);
            indicesPos[q] = rc[0];
            indicesAA[q] = rc[1];
            indicesRot[q] = rc[2];
        }
        
        double energies[][] = stateEMatrix[state].buildReducedMatrix(totUnpruned,indicesPos,indicesAA,indicesRot);
        
        return new ReducedEnergyMatrix(totUnpruned, energies, indicesPos, indicesAA, indicesRot);
    }
    
    
    boolean isSeqFullyDefined(SeqTreeNode seqNode){
        //are we down to a single sequence?
        for(int i=0; i<numTreeLevels; i++){
            if(seqNode.conf[i]>=numAATypes[i])//still a compound option in conf for res i
                return false;
        }
        
        return true;//all single-AA options
    }
    
    
    ConfSpaceSplit nextSplitFixedOrder(SeqTreeNode seqNode){
        //get the next split to use, assuming fixed order (for this tree, and we'll say
        //for ConfTrees too)
        
        if(isSeqFullyDefined(seqNode)){//split a conf tree.  These will be created now
            //to avoid really unbalanced splitting of conf trees,
            //we'll use the somewhat arbitrary heuristic of splitting the earliest firstSplittableLevel 
            //in all the trees (tie broken by using lower-numbered state)
            int firstSplittableLevel = Integer.MAX_VALUE;
            int stateToSplit = -1;
            
            for(int state=0; state<numStates; state++){
                ConfTree stateTree = seqNode.stateTrees[state];
                
                if(stateTree!=null){
                    int stateFirstSplittableLevel = stateTree.getFirstSplittableLevel( stateTree.curExpansion.peek().conf );
                    if( stateFirstSplittableLevel < Math.min(firstSplittableLevel,stateTree.numTreeLevels) ){
                        //tree has a splittable level, and it's the lowest so far
                        stateToSplit = state;
                        firstSplittableLevel = stateFirstSplittableLevel;
                    }
                }
            }
            
            if(stateToSplit==-1)//all trees fully defined!
                return null;
            
            ConfTree treeToSplit = seqNode.stateTrees[stateToSplit];
            return seqNode.stateTrees[stateToSplit].pickSplit( treeToSplit.curExpansion.peek() );
        }
        else {
            for(int i=0; i<numTreeLevels; i++){
                if(seqNode.conf[i]>=numAATypes[i]){//can split this level
                    ArrayList<Integer> opt = new ArrayList<Integer>();//split to all allowed AA types
                    for(int aa=0; aa<numAATypes[i]; aa++)
                        opt.add(aa);
                    return new ConfSpaceSplit(i,opt,null);//null tree indicates sequence split
                }
            }
        }
        
        throw new RuntimeException("ERROR: Sequence wrongly called not fully defined");
    }
    
    
    
    ConfSpaceSplit pickBestSplit(SeqTreeNode seqNode){
        
        if(fixedOrder)
            return nextSplitFixedOrder(seqNode);
        
        throw new RuntimeException("ERROR: dynamic ordering not fully coded yet...");
    }
    
    
    int getStateNum(SeqTreeNode seqNode, ConfTree ct){
        //what state is ct for?
        //It's an error if ct is not one of seqNode's state trees!
        for(int state=0; state<numStates; state++){
            if(seqNode.stateTrees[state]==ct)
                return state;
        }
        
        throw new RuntimeException("ERROR: SeqTree.getStateNum can't find the "
                + "given ConfTree in the given SeqTreeNode's stateTrees!");
    }
    
    
    class GMECLinFunc implements Serializable {
        //this is a function of sequence
        //it is a linear function of the GMEC energies for this sequence in all the states
        //just constTerm + sum_s coeffs_s * GMEC_E_for_state_s

        double[] coeffs;
        double constTerm;
        
        GMECLinFunc(String s){
            //parse from a string listing coefficients in order
            StringTokenizer st = new StringTokenizer(s);
            if(st.countTokens()!=numStates+1){
                throw new RuntimeException("ERROR: SeqTree has "+numStates+" states but GMECLinFunc "
                        + "specified with "+st.countTokens()+" coefficients: "+s);
            }
            
            coeffs = new double[numStates];
            for(int state=0; state<numStates; state++)
                coeffs[state] = Double.valueOf( st.nextToken() );
            
            constTerm = Double.valueOf( st.nextToken() );
        }
    }
    
    
    
    //calculation of node lower bounds
    
    //ok for each obj func and each constr func
    //we can get lower bound using eq 3,4 from seq_sel -- just requires energy matrices
    //and current lowest nodes from each of the trees
    //constraints can be used to prune seqnodes
    //obj func used to decide which to expand though
    //split picking will initially be based on obj func but can handle constr if need be...
    //NEED LIN FUNC OF GMECS AS OBJECT: really just coefficients
    //for constr of form f(seq)<=0, lb(f)>0 for a seqnode allows pruning of node
    //so let's put them in that form
    //since we don't order that way maybe just calc lb and try to prune when expanding nodes...
    
    
    double calcLB(SeqTreeNode seqNode, GMECLinFunc func){
        //get a lower bound for the given GMECLinFunc over the set of sequences allowed in seqNode
        //to avoid repeating code we do this as calcLBSubst with no substitution made
        //(since there is no state -1)
        //return calcLBSubst(seqNode,func,null,-1);
        if(isSeqFullyDefined(seqNode))
            return calcLBConfTrees(seqNode,func);
        else
            return calcLBSubstQuick(seqNode,func,null,-1);
    }
    
    
    double calcLBConfTrees(SeqTreeNode seqNode, GMECLinFunc func){
        //here the sequence is fully defined
        //so we can bound func solely based on lower and upper bounds (depending on func coefficient sign)
        //of the GMECs for each state, which can be derived from the front node of each state's ConfTree
        double ans = func.constTerm;
        for(int state=0; state<numStates; state++){
            if(func.coeffs[state]>0){
                
                if(seqNode.stateTrees[state]==null)//state and sequence impossible
                    return Double.POSITIVE_INFINITY;
                
                ans += func.coeffs[state] * seqNode.stateTrees[state].curExpansion.peek().LB;
            }
            else if(func.coeffs[state]<0){
                
                if(lazyUB){
                    ConfTree curTree = seqNode.stateTrees[state];
                    
                    if(curTree==null){//state and sequence impossible
                        //bound would be +infinity
                        ans = Double.NEGATIVE_INFINITY;
                        continue;
                    }
                    
                    //make sure stateUB is updated, at least based on the current best node in this state's tree
                    ConfTreeNode curNode = curTree.curExpansion.peek();
                    
                    if(Double.isNaN(curNode.UB)){//haven't calculated UB, so calculate and update stateUB
                        while(!curTree.updateUB(curNode, curNode.UBConf)){
                            //if UB not calc'd yet, curNode.UBConf is from curNode's parent
                            //if updateUB fails then the node has an inevitable clash...remove it from the tree
                            curTree.curExpansion.poll();
                            if(curTree.curExpansion.isEmpty()){//no more nodes, so no available states!
                                //bound would be +infinity
                                seqNode.stateUB[state] = Double.POSITIVE_INFINITY;
                                ans = Double.NEGATIVE_INFINITY;
                                continue;
                            }
                            
                            curNode = curTree.curExpansion.peek();
                            if(!Double.isNaN(curNode.UB))
                                break;
                        }
                    }
                    
                    seqNode.stateUB[state] = Math.min(seqNode.stateUB[state],curNode.UB);    
                }
                
                ans += func.coeffs[state] * seqNode.stateUB[state];
            }
            
            //shell-shell energies may differ between states!
            ans += func.coeffs[state]*stateEMatrix[state].getShellShellE();
        }
        
        return ans;
    }
    
        
    //quick non-LP bound
    double calcLBSubstQuick(SeqTreeNode seqNode, GMECLinFunc func, ConfTreeNode altBestConfNode, int substState){
        
        double ans = func.constTerm;
        //first terms for mutable residues
        for(int i=0; i<numTreeLevels; i++){
            double resE = Double.POSITIVE_INFINITY;
            
            for( int curAA : AATypeOptions.get(i).get(seqNode.conf[i]) ){
                double AAE = 0;
                
                //get contributions to residue energy for this AA type from all states
                for(int state=0; state<numStates; state++){
                    if(func.coeffs[state]!=0){
                        
                        int stateResNum = mutable2StateResNums.get(state).get(i);//residue number i converted to this state's flexible residue numbering
                        boolean minForState = (func.coeffs[state]>0);//minimizing (instead of maximizing) energy for this state

                        double stateAAE = Double.POSITIVE_INFINITY;

                        ArrayList<Integer> rotList = stateRCsAtAA.get(state).get(stateResNum).get(curAA);
                        PrunedRotamers<Boolean> prunedRot = seqNode.statePrunedRot[state].eliminatedRotAtPos;
                        boolean[][][][][][] splitFlags = seqNode.statePrunedPairs[state].splitFlags;
                        PairwiseEnergyMatrix eMatrix = stateEMatrix[state];

                        if(!minForState){
                            stateAAE = Double.NEGATIVE_INFINITY;
                            rotList = negStateRCsAtAA.get(state).get(stateResNum).get(curAA);
                            prunedRot = seqNode.negStatePrunedRot[state].eliminatedRotAtPos;
                            splitFlags = seqNode.negStatePrunedPairs[state].splitFlags;
                            eMatrix = negStateEMatrix[state];
                        }

                        for(int rot : rotList){
                            //make sure rot isn't pruned
                            if(!prunedRot.get(stateResNum, curAA, rot)){

                                double rotE = eMatrix.getIntraAndShellE(stateResNum,curAA,rot);

                                for(int res2=0; res2<stateNumRes[state]; res2++){//all non-mut; seq only if < this one
                                    if( (!mutable2StateResNums.get(state).contains(res2)) || res2<stateResNum){

                                        double bestInteraction = Double.POSITIVE_INFINITY;
                                        ArrayList<int[]> rotList2 = stateRCs.get(state).get(res2);

                                        if(!minForState){
                                            bestInteraction = Double.NEGATIVE_INFINITY;
                                            rotList2 = negStateRCs.get(state).get(res2);
                                        }

                                        for(int rot2[] : rotList2){
                                            if(!prunedRot.get(res2, rot2[0], rot2[1])){
                                                if(!splitFlags[stateResNum][curAA][rot][res2][rot2[0]][rot2[1]]){

                                                    double pairwiseE = eMatrix.getPairwiseE(stateResNum,curAA,rot,res2,rot2[0],rot2[1]);

                                                    if(minForState)
                                                        bestInteraction = Math.min(bestInteraction,pairwiseE);
                                                    else
                                                        bestInteraction = Math.max(bestInteraction,pairwiseE);
                                                }
                                            }
                                        }

                                        rotE += bestInteraction;
                                    }
                                }

                                if(minForState)
                                    stateAAE = Math.min(stateAAE,rotE);
                                else
                                    stateAAE = Math.max(stateAAE,rotE);
                            }
                        }

                        if(Double.isInfinite(stateAAE)){
                            //this will occur if the state is impossible (all confs pruned)
                            if(func.coeffs[state]>0)
                                AAE = Double.POSITIVE_INFINITY;
                            else if(func.coeffs[state]<0 && AAE!=Double.POSITIVE_INFINITY)
                                //if a "positive-design" (coeff>0) state is impossible, return +inf overall
                                AAE = Double.NEGATIVE_INFINITY;
                            //else AAE unchanged, since func doesn't involve this state
                        }
                        else
                            AAE += func.coeffs[state]*stateAAE;
                    }
                }
                
                resE = Math.min(resE,AAE);
            }
            
            ans += resE;
        }
        
        //now we bound the energy for the other residues for each of the states
        //(internal energy for that set of residues)
        for(int state=0; state<numStates; state++){
            
            if(func.coeffs[state]!=0){
                ans += func.coeffs[state]*stateEMatrix[state].getShellShellE();//shell-shell E may vary btw states...
                //though should be the same for positive, negative versions of a state

                double nonMutBound = boundStateNonMutE(state,seqNode,func.coeffs[state]>0);
                //handle this like AAE above
                if(Double.isInfinite(nonMutBound)){
                    //this will occur if the state is impossible (all confs pruned)
                    if(func.coeffs[state]>0)
                        return Double.POSITIVE_INFINITY;
                    else if(func.coeffs[state]<0)
                        ans = Double.NEGATIVE_INFINITY;
                        //if a "positive-design" (coeff>0) state is impossible, return +inf overall still
                    //else ans unchanged, since func doesn't involve this state
                }
                else
                    ans += func.coeffs[state]*nonMutBound;
            }
        }
        
        return ans;
    }
    
    
    
    double boundStateNonMutE(int state, SeqTreeNode seqNode, boolean minForState){
        //get a quick lower or upper bound (as indicated) for the energy of the given state's
        //non-mutable residues (their intra+shell energies + pairwise between them)
        //use pruning information from seqNode
        double ans = 0;
        
        for(int res=0; res<stateNumRes[state]; res++){
            if((!mutable2StateResNums.get(state).contains(res))){   
                
                double resE = Double.POSITIVE_INFINITY;
                
                ArrayList<int[]> rotList = stateRCs.get(state).get(res);
                PrunedRotamers<Boolean> prunedRot = seqNode.statePrunedRot[state].eliminatedRotAtPos;
                boolean[][][][][][] splitFlags = seqNode.statePrunedPairs[state].splitFlags;
                PairwiseEnergyMatrix eMatrix = stateEMatrix[state];
                
                if(!minForState){
                    resE = Double.NEGATIVE_INFINITY;
                    rotList = negStateRCs.get(state).get(res);
                    prunedRot = seqNode.negStatePrunedRot[state].eliminatedRotAtPos;
                    splitFlags = seqNode.negStatePrunedPairs[state].splitFlags;
                    eMatrix = negStateEMatrix[state];
                }
                
                for(int rot[] : rotList){
                    //make sure rot isn't pruned
                    if(!prunedRot.get(res, rot[0], rot[1])){

                        double rotE = eMatrix.getIntraAndShellE(res, rot[0], rot[1]);

                        for(int res2=0; res2<stateNumRes[state]; res2++){//all non-mut; seq only if < this one
                            if( (!mutable2StateResNums.get(state).contains(res2)) && res2<res){

                                double bestInteraction = Double.POSITIVE_INFINITY;
                                ArrayList<int[]> rotList2 = stateRCs.get(state).get(res2);
                                if(!minForState){
                                    bestInteraction = Double.NEGATIVE_INFINITY;
                                    rotList2 = negStateRCs.get(state).get(res2);
                                }
                                
                                for(int rot2[] : rotList2){
                                    if(!prunedRot.get(res2, rot2[0], rot2[1])){
                                        if(!splitFlags[res][rot[0]][rot[1]][res2][rot2[0]][rot2[1]]){

                                            double pairwiseE = eMatrix.getPairwiseE(res,rot[0],rot[1],res2,rot2[0],rot2[1]);

                                            if(minForState)
                                                bestInteraction = Math.min(bestInteraction,pairwiseE);
                                            else
                                                bestInteraction = Math.max(bestInteraction,pairwiseE);
                                        }
                                    }
                                }

                                rotE += bestInteraction;
                            }
                        }

                        if(minForState)
                            resE = Math.min(resE,rotE);
                        else
                            resE = Math.max(resE,rotE);
                    }
                }
                
                ans += resE;
            }
        }
        
        return ans;
    }
    
    
    
    //regenerate DEEGoldsteinPairs from NewPrunedPairs (array for residue pairs) for given state
    //the single pruned rots and triples to use are provided
    DEEGoldsteinPairs uncompressPrunedPairs(NewPrunedPairs[][] npp,
            PrunedRotamers<Boolean> prunedRot, int state,
            boolean[][][][][][][][][] tripFlags){
        
        
        DEEGoldsteinPairs template = npp[0][1].getFullAncestor();
        
        //regenerate splitFlags
        //initialize (all false)
        boolean[][][][][][] splitFlags = template.pairwiseMinEnergyMatrix.initializePairwiseBooleanMatrix();

        
        for(int res1=0; res1<stateNumRes[state]; res1++){
            for(int res2=0; res2<stateNumRes[state]; res2++){
                if(res1!=res2){
                    npp[res1][res2].fillInSplitFlags(splitFlags,res1,res2);
                }
            }
        }
               
        
        return new DEEGoldsteinPairs( template.pairwiseMinEnergyMatrix, template.pairwiseMaxEnergyMatrix,
                        template.numMutable, template.strandMut, template.curEw-template.Ival, 
                        template.strandRot, prunedRot, template.resInPair, template.doMinimize, splitFlags, 
                        template.useFlags, template.magicBullet, template.distrDEE, template.minimizeBB, 
                        template.maxScale!=1, template.maxScale, template.mutRes2Strand, template.mutRes2MutIndex, 
                        template.typeDependent, template.doIMinDEE,
                        template.Ival, tripFlags, template.doPerturbations, false );
    }
    
    
    void compressPrunedPairs(SeqTreeNode seqNode, SeqTreeNode parent){
        //compress the pruned pairs to NewPrunedPairs
        //if parent is null, then there's already a NewPrunedPairs here and we can overwrite
        //using the same parent NewPrunedPairs
        //otherwise we are talking about new pairs pruned since parent
                
        NewPrunedPairs newNPP[][][] = new NewPrunedPairs[numStates][][];
        NewPrunedPairs newNegNPP[][][] = new NewPrunedPairs[numStates][][];
        
        for(int state=0; state<numStates; state++){
            
            NewPrunedPairs[][] parentNPP = null;
            DEEGoldsteinPairs parentDGP = null;
        
            if(parent==null){
                //updating node
                //so it should already have npp, just need to update with any new pruning...
                parentNPP = seqNode.npp[state];
            }
            else if(parent.npp==null)
                parentDGP = parent.statePrunedPairs[state];
            else
                parentNPP = parent.npp[state];
                
            
            newNPP[state] = new NewPrunedPairs[stateNumRes[state]][stateNumRes[state]];
            
            for(int res1=0; res1<stateNumRes[state]; res1++){
                for(int res2=0; res2<stateNumRes[state]; res2++){
                    if(res1!=res2){
                        
                        if(parentDGP!=null){
                            newNPP[state][res1][res2] = 
                                    new NewPrunedPairs(stateRCs.get(state).get(res1),
                                    stateRCs.get(state).get(res2),res1,res2,
                                    seqNode.statePrunedPairs[state].splitFlags,parentDGP,null);
                        }
                        else {
                            newNPP[state][res1][res2] = 
                                    new NewPrunedPairs(stateRCs.get(state).get(res1),
                                    stateRCs.get(state).get(res2),res1,res2,
                                    seqNode.statePrunedPairs[state].splitFlags,null,
                                    parentNPP[res1][res2]);
                        }
                    }
                }
            }
            
            if(seqNode.statePrunedPairs[state] == seqNode.negStatePrunedPairs[state])
                newNegNPP[state] = newNPP[state];
            else {//have to compress negative state separately
                parentNPP = null;
                parentDGP = null;

                if(parent==null){
                    //updating node
                    //so it should already have npp, just need to update with any new pruning...
                    parentNPP = seqNode.negNPP[state];
                }
                else if(parent.negNPP==null)
                    parentDGP = parent.negStatePrunedPairs[state];
                else
                    parentNPP = parent.negNPP[state];


                newNegNPP[state] = new NewPrunedPairs[stateNumRes[state]][stateNumRes[state]];

                for(int res1=0; res1<stateNumRes[state]; res1++){
                    for(int res2=0; res2<stateNumRes[state]; res2++){
                        if(res1!=res2){

                            if(parentDGP!=null){
                                newNegNPP[state][res1][res2] = 
                                        new NewPrunedPairs(negStateRCs.get(state).get(res1),
                                        negStateRCs.get(state).get(res2),res1,res2,
                                        seqNode.negStatePrunedPairs[state].splitFlags,parentDGP,null);
                            }
                            else {
                                newNegNPP[state][res1][res2] = 
                                        new NewPrunedPairs(negStateRCs.get(state).get(res1),
                                        negStateRCs.get(state).get(res2),res1,res2,
                                        seqNode.negStatePrunedPairs[state].splitFlags,null,
                                        parentNPP[res1][res2]);
                            }
                        }
                    }
                }
            }
        }
        
        
        seqNode.npp = newNPP;
        seqNode.negNPP = newNegNPP;

        //and now we get rid of the uncompressed pairs
        seqNode.statePrunedPairs = null;
        seqNode.negStatePrunedPairs = null;
        //the split flags are also in the rot-pruning objects
        for(int state=0; state<numStates; state++){
            seqNode.statePrunedRot[state].splitFlags = null;
            seqNode.negStatePrunedRot[state].splitFlags = null;
        }
    }
    
    

    
}
