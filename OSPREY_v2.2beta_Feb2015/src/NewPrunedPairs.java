
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;

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
//	NewPrunedPairs.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class NewPrunedPairs implements Serializable {
    //compactly represents newly pruned pairs at some SeqTreeNode
    //at some residue pair
    
    HashSet<Integer> newPairs;
    int maxIndex1;
    
    //these new pairs are added onto either a NewPrunedPairs object
    //(to indicate recursive pruning on many levels)
    //or onto a full DEEGoldsteinPairs object (so one of these should be null)
    DEEGoldsteinPairs parentDGP = null;
    NewPrunedPairs parentNPP = null;
    
    //these are stateRCs lists for our residues of interests
    //probably shared with SeqTree
    ArrayList<int[]> RCList1, RCList2;
    
    
    public NewPrunedPairs(ArrayList<int[]> RCList1, ArrayList<int[]> RCList2, 
            int res1, int res2,
            boolean[][][][][][] splitFlags, DEEGoldsteinPairs parentDGP, NewPrunedPairs parentNPP){
        //(AA,rot) lists for the RCs at the res involved, which res they are, 
        //and the splitFlags array to compress
        //also the higher-level stuff to compare to
        
        this.RCList1 = RCList1;
        this.RCList2 = RCList2;
        this.parentDGP = parentDGP;
        this.parentNPP = parentNPP;
        
        if( (parentDGP==null) == (parentNPP==null) ){
            throw new RuntimeException("ERROR: To set up NewPrunedPairs"
                    + " parent must be exactly one of DEEGoldsteinPairs or NewPrunedPairs");
        }

        maxIndex1 = RCList1.size();
        
        newPairs = new HashSet<>();
        
        for(int r1=0; r1<RCList1.size(); r1++){
            for(int r2=0; r2<RCList2.size(); r2++){
                int rot1[] = RCList1.get(r1);
                int rot2[] = RCList2.get(r2);
                if(splitFlags[res1][rot1[0]][rot1[1]][res2][rot2[0]][rot2[1]]){//If pair pruned in splitFlags...
                    if(!isPrunedByAncestors(res1,r1,res2,r2)){//but not in ancestors...
                        newPairs.add(singleIndex(r1,r2));
                    }
                }
            }
        }
    }
    
    
    void fillInSplitFlags(boolean[][][][][][] splitFlags, int res1, int res2){
        //fill in the split flags for the specified residues based on this
        for(int r1=0; r1<RCList1.size(); r1++){
            for(int r2=0; r2<RCList2.size(); r2++){
                int rot1[] = RCList1.get(r1);
                int rot2[] = RCList2.get(r2);
                
                if(isPrunedByAncestors(res1,r1,res2,r2))
                    splitFlags[res1][rot1[0]][rot1[1]][res2][rot2[0]][rot2[1]] = true;
                else if(isPrunedHere(r1,r2))
                    splitFlags[res1][rot1[0]][rot1[1]][res2][rot2[0]][rot2[1]] = true;
            }
        }
    }
    
    
    
    DEEGoldsteinPairs getFullAncestor(){
        //somewhere along the chain of descent there should be an uncompressed DEEGoldsteinPairs
        //find it
        if(parentDGP==null)
            return parentNPP.getFullAncestor();
        else
            return parentDGP;
    }
    
    int singleIndex(int r1, int r2){
        //convert rotamer indices to single index
        //we chose maxIndex1 so this mapping is 1-to-1
        
        //ensure no overflow
        if(Integer.MAX_VALUE/(r2+1)<=maxIndex1)
            throw new RuntimeException("ERROR: Overflow in newPrunedPairs.singleIndex.  max int: "
                    +Integer.MAX_VALUE+" maxIndex1: "+maxIndex1+" r2: "+r2);
        
        return r1 + r2*maxIndex1;
    }
    
    
    boolean isPrunedHere(int r1, int r2){
        return newPairs.contains(singleIndex(r1,r2));
    }
                
        
    
    
    boolean isPrunedByAncestors(int res1, int r1, int res2, int r2){
        if(parentDGP!=null){
            //ancestor is just the parent pruned pairs...look up there
            
            int rot1[] = RCList1.get(r1);
            int rot2[] = RCList2.get(r2);
            if(parentDGP.splitFlags[res1][rot1[0]][rot1[1]][res2][rot2[0]][rot2[1]])
                return true;
            else
                return false;
        }
        else {
            //see if it's pruned in the parent
            if(parentNPP.isPrunedHere(r1, r2))
                return true;
            else//otherwise see if it's pruned in any older ancestors
                return parentNPP.isPrunedByAncestors(res1, r1, res2, r2);
        }
    }
    
    
}
