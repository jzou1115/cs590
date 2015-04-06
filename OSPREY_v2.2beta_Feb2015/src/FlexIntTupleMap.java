
import java.util.ArrayList;

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
//	FlexIntTupleMap.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class FlexIntTupleMap<T> {
    //like IntTupleMap but allowing different tuple size
    //we achieve this by having an IntTupleMap for each tuple size
    IntTupleMap[] maps;//indexed by tuple size
    
    public FlexIntTupleMap(int... top){
        maps = new IntTupleMap[top.length+1];
        for(int a=0; a<=top.length; a++){
            int[] redTop = new int[a];
            System.arraycopy(top, 0, redTop, 0, a);
            maps[a] = new IntTupleMap(redTop);
        }
    }
    
    
    void set(T t, int... indices){
        maps[indices.length].set(t,indices);
    }
    
    void remove(int... indices){
        maps[indices.length].remove(indices);
    }
    
    T get(int... indices){
        return (T)maps[indices.length].get(indices);
    }
    
    ArrayList<int[]> getKeys(){
        ArrayList<int[]> ans = new ArrayList<>();
        for(IntTupleMap map : maps){
            ans.addAll(map.getKeys());
        }
        return ans;
    }
    
    
    
}
