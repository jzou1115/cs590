
import java.io.Serializable;

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
//	ConfTreeNode.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

//A node to go in the expansion for a ConfTree
//represents a subset of conformational space
//sorted in ConfTree's priority queue using fScore
class ConfTreeNode implements Comparable, Serializable {
    
    int conf[];//conformational space expressed as the current compound option (defined in confTree) at each level
    //-1 means all options possible at the given level
    
    //lower and upper bounds on the minimum energy for the conformational space represented by this node:
    double LB = Double.NEGATIVE_INFINITY;//current lower bound.  Initialized with ConfTree.quickFScore
    //can be tightened lazily later, e.g., using EPIC to represent continuous terms
    
    double UB = Double.POSITIVE_INFINITY;//upper bound
    int UBConf[] = null;//can have an upper bound on GMEC energy for this node's conf space
    //(and thus on the overall GMEC energy)
    //that is the energy of the conf denoted by UBConf (consisting entirely of simple options) 
    
    //comparison is by lower bound, for use in A*
    //COULD SOMETIMES WANT BY UB TOO
    @Override
    public int compareTo (Object otherObject) {
        
        ConfTreeNode other = (ConfTreeNode)otherObject;
        int comparison = Double.valueOf(LB).compareTo(other.LB);
        
        if(comparison==0){
            //need to break tie in LB's
            //let's go through the conformation lexically
            //we can assume confs are of same length (number of levels for the ConfTree)
            for(int level=0; level<conf.length; level++){
                comparison = Integer.valueOf(conf[level]).compareTo(other.conf[level]);
                if(comparison!=0)
                    return comparison;
            }
            
            return 0;//same conf--thus should be equal
        }
        else
            return comparison;
    }
        
    
}
