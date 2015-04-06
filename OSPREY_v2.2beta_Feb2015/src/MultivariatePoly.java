
import cern.colt.matrix.DoubleMatrix1D;
import java.util.ArrayList;
import java.util.Arrays;

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
//	MultivariatePoly.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class MultivariatePoly {
    //A multivariate polynomial
    //pretty basic for now...can be optimized if it becomes a bottleneck
    
    FlexIntTupleMap<Double> terms;//map variable #s in descending order to their coefficient
    int numDOFs;
    int order;//maximum total order of polynomial
    
    
    
    
    public MultivariatePoly(int numDOFs, int order){
        //start with 0 polynomial
        init(numDOFs,order);
    }
      
    
    void init(int numDOFs, int order){
        int[] termBounds = new int[order];
        Arrays.fill(termBounds, numDOFs);
        this.numDOFs = numDOFs;
        this.order = order;
        terms = new FlexIntTupleMap<>(termBounds);
    }
    
    
    static MultivariatePoly one(int numDOFs){//constant polynomial equal to 1
        MultivariatePoly ans = new MultivariatePoly(numDOFs,0);
        ans.terms.set(1.);
        return ans;
    }
    
    
    //conversion to and from SeriesFitter coeff array
    //adapted from functions in SeriesFitter (evalSeries, etc.)
    public MultivariatePoly(double coeffs[], int nd, int order, boolean includeConst, int PCOrder, boolean isPC[]){
        init(nd,order);
        
        //we support orders 2-6
        //and PCOrder may range up to 6 (but has no effect if <=order)
        if(order<2||order>6||PCOrder>6){
            throw new RuntimeException("ERROR: MultivariatePoly constructor does not support order "
                    +order+" and/or PCOrder "+PCOrder);
        }
        
        
        int count = 0;

        if(includeConst){
            terms.set(coeffs[0]);
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            terms.set(coeffs[count],dof);
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            for(int dof2=0; dof2<dof; dof2++){
                terms.set(coeffs[count],dof,dof2);
                count++;
            }
            terms.set(coeffs[count],dof,dof);
            count++;
        }


        if(order>=3){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){

                        terms.set(coeffs[count],dof,dof2,dof3);
                        count++;
                    }
                }
            }
        }
        else if(PCOrder>=3){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    terms.set(coeffs[count],dof,dof2,dof3);
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            terms.set(coeffs[count],dof,dof2,dof3,dof4);
                            count++;
                        }
                    }
                }
            }
        }
        else if(PCOrder>=4){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            terms.set(coeffs[count],dof,dof2,dof3,dof4);
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        if(order>=5){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    terms.set(coeffs[count],dof,dof2,dof3,dof4,dof5);
                                                    count++;
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=5){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    terms.set(coeffs[count],dof,dof2,dof3,dof4,dof5);
                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=6){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                            terms.set(coeffs[count],dof,dof2,dof3,dof4,dof5,dof6);
                                                            count++;
                                                    }
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=6){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        if(isPC[dof6]){

                                                            terms.set(coeffs[count],dof,dof2,dof3,dof4,dof5,dof6);
                                                            count++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
    }
    
    
    
    double[] toSeriesDoubleArray(boolean includeConst, int PCOrder, boolean isPC[]){

        double coeffs[] = new double[SeriesFitter.getNumParams(numDOFs, includeConst, order)];
        int nd = numDOFs;
        
        //we support orders 2-6
        //and PCOrder may range up to 6 (but has no effect if <=order)
        if(order<2||order>6||PCOrder>order){
            throw new RuntimeException("ERROR: MultivariatePoly.toSeriesDoubleArray does not support order "
                    +order+" and/or PC");
        }
        
        int count = 0;

        if(includeConst){
            coeffs[0] = terms.get();
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            coeffs[count] = terms.get(dof);
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            for(int dof2=0; dof2<dof; dof2++){
                coeffs[count] = terms.get(dof,dof2);
                count++;
            }
            coeffs[count] = terms.get(dof,dof);
            count++;
        }


        if(order>=3){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){

                        coeffs[count] = terms.get(dof,dof2,dof3);
                        count++;
                    }
                }
            }
        }
        else if(PCOrder>=3){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    coeffs[count] = terms.get(dof,dof2,dof3);
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            coeffs[count] = terms.get(dof,dof2,dof3,dof4);
                            count++;
                        }
                    }
                }
            }
        }
        else if(PCOrder>=4){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            coeffs[count] = terms.get(dof,dof2,dof3,dof4);
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        if(order>=5){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    coeffs[count] = terms.get(dof,dof2,dof3,dof4,dof5);
                                                    count++;
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=5){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    coeffs[count] = terms.get(dof,dof2,dof3,dof4,dof5);
                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=6){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                            coeffs[count] = terms.get(dof,dof2,dof3,dof4,dof5,dof6);
                                                            count++;
                                                    }
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=6){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        if(isPC[dof6]){

                                                            coeffs[count] = terms.get(dof,dof2,dof3,dof4,dof5,dof6);
                                                            count++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return coeffs;
    }
    
    
    //arithmetic on polynomials
    void add(MultivariatePoly pol2){
        //add pol2 to this
        for(int[] term2 : pol2.terms.getKeys()){
            Double coeff1 = terms.get(term2);
            if(coeff1==null)
                terms.set((double)pol2.terms.get(term2),term2);
            else
                terms.set(pol2.terms.get(term2)+coeff1,term2);
        }
    }
    
    
    MultivariatePoly times(MultivariatePoly pol2){
        //return this*pol2
        MultivariatePoly ans = new MultivariatePoly(numDOFs,order+pol2.order);
        for(int[] term1 : terms.getKeys()){
            for(int[] term2 : pol2.terms.getKeys()){
                
                int[] prodTerm = new int[term1.length+term2.length];//merge of term1 and term2
                int count1=0, count2=0;
                for(int a=0; a<prodTerm.length; a++){
                    if(count1==term1.length){
                        prodTerm[a] = term2[count2];
                        count2++;
                    }
                    else if(count2==term2.length){
                        prodTerm[a] = term1[count1];
                        count1++;
                    }
                    else if(term1[count1]>term2[count2]){
                        prodTerm[a] = term1[count1];
                        count1++;
                    }
                    else {
                        prodTerm[a] = term2[count2];
                        count2++;
                    }
                }
                    
                Double coeff1 = ans.terms.get(prodTerm);
                double prodCoeff = terms.get(term1)*pol2.terms.get(term2);
                
                
                if(coeff1==null)
                    ans.terms.set(prodCoeff,prodTerm);
                else
                    ans.terms.set(prodCoeff+coeff1,prodTerm);
            }
        }
        
        return ans;
    }
    
    
    
    MultivariatePoly shift(DoubleMatrix1D c){
        //if this is f(x), return f(x+c)
        MultivariatePoly ans = new MultivariatePoly(numDOFs,order);
        
        MultivariatePoly newVars[] = new MultivariatePoly[numDOFs];
        //x_0+c+0, x_1+c_1,...
        for(int dof=0; dof<numDOFs; dof++){
            newVars[dof] = new MultivariatePoly(numDOFs,1);
            newVars[dof].terms.set(c.get(dof));
            newVars[dof].terms.set(1.,dof);
        }

        for(int[] term1 : terms.getKeys()){
            //expand each monomial of this
            MultivariatePoly expansion = new MultivariatePoly(numDOFs,0);
            expansion.setConstant(terms.get(term1));
            for(int dof : term1)
                expansion = expansion.times(newVars[dof]);
            ans.add(expansion);
        }
        
        return ans;
    }
    
    
    double getConstant(){
        //return constant term
        return terms.get();//get [] that is
    }
    
    void removeConstant(){
        //remove constant terms
        terms.remove();
    }
    
    void setConstant(double a){
        //set constant term
        terms.set(a);
    }
    
    void print(){
        ArrayList<int[]> keys = terms.getKeys();
        for(int k=0; k<keys.size(); k++){
            if(k!=0)
                System.out.print("+");
            System.out.print(terms.get(keys.get(k)));
            for(int dof : keys.get(k))
                System.out.print("*x"+dof);
        }
        System.out.println();
    }
}