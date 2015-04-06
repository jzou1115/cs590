/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

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
//	SeqTreeNode.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
class SeqTreeNode extends ConfTreeNode {
    
    //we have conf (here meaning allowed AA-type sets at each res), UB and LB (here meaning bounds on GMEC 
    //for sequence space), inherited from ConfTreeNode.  And we compare by LB in the same way
        
    ConfTree[] stateTrees;//conformational trees for each of the states
    //if a state doesn't need to be LB'd and the sequence is not fully defined, we don't need this yet...
    
    
    //for states where we need to upper-bound conf, also need to keep updating pruned rots/tuples
    //these classes can update these as needed and also keep track of them
    DEEGoldstein statePrunedRot[];
    
    DEEGoldsteinPairs statePrunedPairs[] = null;
    NewPrunedPairs[][][] npp = null;//compressed pairs: only list new pruned pairs at this level
    //SeqTree can perform (de)compression
    
    //DEEGoldsteinTriples statePrunedTriples[];//maybe later
    
    
    double stateUB[];//if this node's sequence is fully defined,
    //stateUB will upper-bound the GMEC for each state
    
    //may also want pruning for negative versions of states (e.g. discrete instead of continuous)
    //can be set equal to the normal versions
    DEEGoldstein negStatePrunedRot[];
    DEEGoldsteinPairs negStatePrunedPairs[];
    NewPrunedPairs[][][] negNPP = null;
    
    
    void setupPruning(RotamerSearch rs, boolean neg, int state, int[][][] stateStrandMut){
        //set up pruning structures
        //uncompressed for now
        //called from SeqTree constructor
        
        
        //this stuff is suitable for either rigid, or traditional minDEE
        //if need to change can use Multistate.cfg
        double Ew = 0;
        boolean minimizeBB = false;
        boolean imindee = false;
        double Ival = 0;
        boolean scaleInt = false;
        double maxScale = 0;
        
        //minimization in SeqTree is by EPIC so check for cetm to see if rs minimizes
        boolean doMin = (rs.cetm!=null);


        DEEGoldstein prunedRot = new DEEGoldstein(rs.getMinMatrix(), 
                rs.getMaxMatrix(), rs.numberMutable, stateStrandMut[state], Ew, 
                rs.strandRot, rs.eliminatedRotAtRes, doMin,
                rs.indIntMinDEE, rs.pairIntMinDEE, rs.splitFlags, true, minimizeBB, 
                rs.mutRes2Strand, rs.mutRes2StrandMutIndex, true, imindee, 
                Ival, rs.doPerturbations,false);

        //not doing triples for now
        boolean[][][][][][][][][] tripleFlags = null;

        DEEGoldsteinPairs prunedPairs = new DEEGoldsteinPairs(rs.getMinMatrix(), 
                rs.getMaxMatrix(), rs.numberMutable, stateStrandMut[state], Ew, 
                rs.strandRot, rs.eliminatedRotAtRes, null, doMin, 
                rs.splitFlags, true, false/*magicBullet*/, false, minimizeBB, 
                scaleInt, maxScale, rs.mutRes2Strand, rs.mutRes2StrandMutIndex, true, imindee,
                Ival, tripleFlags, rs.doPerturbations, false );
        
        if(neg) {
            negStatePrunedRot[state] = prunedRot;
            negStatePrunedPairs[state] = prunedPairs;
        }
        else {
            statePrunedRot[state] = prunedRot;
            statePrunedPairs[state] = prunedPairs;
        }
    }
}