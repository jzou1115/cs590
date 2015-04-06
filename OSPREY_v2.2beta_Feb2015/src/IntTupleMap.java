/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
//This is taken from gnn
import java.util.ArrayList;
import java.util.TreeMap;

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
//	IntTupleMap.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class IntTupleMap<T> {
    //map from tuple of nonnegative integers to T
    //the integers are strictly bounded above by indexTop
    
    int indexTop[];
    int numIndices;//size of tuple
    TreeMap<Integer,T> map = new TreeMap<>();
    
    public IntTupleMap(int... top){
        indexTop = top;
        numIndices = indexTop.length;
        
        //check that we can fit everything into a single integer
        int multiplier = 1;
        for(int a=0; a<numIndices; a++){
            if( Integer.MAX_VALUE/multiplier < indexTop[a] )
                throw new RuntimeException("ERROR: Integer overflow for IntTupleMap!");
            else if( indexTop[a]<=0 )
                throw new RuntimeException("ERROR: indexTop must be all positive, not "+indexTop[a]);
            else
                multiplier *= indexTop[a];
        }
    }
    
    
    void set(T t, int... indices){
        map.put(toSingleIndex(indices),t);
    }
    
    void remove(int... indices){
        map.remove(toSingleIndex(indices));
    }
    
    T get(int... indices){
        int singleIndex = toSingleIndex(indices);
        if(map.containsKey(singleIndex))
            return map.get(singleIndex);
        else
            return null;
    }
    
    ArrayList<int[]> getKeys(){
        ArrayList<int[]> ans = new ArrayList<>();
        for(int singleIndex : map.keySet()){
            ans.add(toTuple(singleIndex));
        }
        return ans;
    }
    
    
    ArrayList<T> getList(int... indices){
        //like get, except some of the indices may be -1
        //which we treat as a wildcard, returning anything that matches
        
        ArrayList<T> ans = new ArrayList<>();
        
        for(int a=0; a<numIndices; a++){//look for -1's
            if(indices[a]==-1){
                for(int ind=0; ind<indexTop[a]; ind++){
                    int[] indices2 = indices.clone();
                    indices2[a] = ind;
                    ans.addAll(getList(indices2));
                }
                return ans;
            }
        }
        
        //if we get here, there are no -1's in indices
        T t = get(indices);
        if(t!=null)
            ans.add(t);
        
        return ans;
    }
    
    
    int toSingleIndex(int[] tup){//wrap up tup into single index so it can be mapped
        //input checking
        if(tup.length!=numIndices)
            throw new RuntimeException("Error: wrong number of indices in IntTupleMap!");
        
        int ans = 0;
        int multiplier = 1;
        
        for(int a=0; a<numIndices; a++){
            if(tup[a]<0||tup[a]>=indexTop[a])
                throw new RuntimeException("Error: Index "+a+"="+tup[a]+" out of bounds for IntTupleMap!");
            
            ans += multiplier*tup[a];
            multiplier *= indexTop[a];
        }
        
        return ans;
    }
    
    //reverse operation
    int[] toTuple(int singleIndex){
        if(numIndices==0)
            return new int[] {};//just one possible tuple
        
        int[] ans = new int[numIndices];
        for(int a=0; a<numIndices-1; a++){
            ans[a] = singleIndex % indexTop[a];
            singleIndex /= indexTop[a];
        }
        ans[numIndices-1] = singleIndex;
        
        return ans;
    }
    
    
    
}
