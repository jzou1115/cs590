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
//	PartialMinConstrainer.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////


import cern.colt.matrix.DoubleMatrix1D;
import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.linear.LinearConstraint;
import org.apache.commons.math3.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optimization.linear.Relationship;
import org.apache.commons.math3.optimization.linear.SimplexSolver;
import org.apache.commons.math3.optimization.linear.UnboundedSolutionException;


public class PartialMinConstrainer {
    //Get tighter lower bounds on a conformation's energy by constraining partial minima
    
    //maybe we can get tighter constraints on some sums of terms by DEE-type methods too
    //(i.e., show anything better impossible by comparison to competitor(s))
    
    PairwiseEnergyMatrix eMatrix;
    //PairwiseEnergyMatrix pairE;//the pairwise entries in this are minimized full-pair energies

    CETMatrix cetm;//if null then just do discrete portion

    
    
    
    static boolean scaleResDom = true;//scale down other intra+shell terms in res-dom by (numRes-1)
    //so that the res-dominant energies add up to exactly twice the total energy
    //this way we don't get loosening of bounds because of low energy being transferred to pairwise terms
    
    
    boolean[][][][][][] splitFlags;
    boolean[][][][][][][][][] tripleFlags;
    //we expect to be given only unpruned (AAind,rot) but these may have pruned pairs or triples
    
    
    //whether to use larger tuples for dominant undecided-tuple constr
    boolean do2Res = true, do3Res= true;
    
    
    //we'll include single and pair constr regardless (easy)
    //but then we can enhance with either undecided-tuple constr (decidedBased=true)
    //which is easy to port to continuous
    //else go with star and dominant
    boolean decidedBased = true;

    boolean doTupComb = true;
    //this applies constraints for
    //each combination of energies involving up to a certain number of undecided residues
    //(decided always all included)
    //CURRENTLY IF APPLIED EXCLUDES ALL OTHERS
    //LATER TRY TO COMBINE WITH OTHER CONSTRAINTS MAYBE?
    
    //constrain all combinations of energy terms within a tuple of interest
    //for this purpose we will group all undecided together
    
    
    
    boolean contOnly = false;//can be used to bound only the continuous contribution
    
    //LETS SET UP INEQ LIST WITH CONSTRUCTOR
    //ADD INEQS W/ DEDICATED FUNCTION
    ArrayList<LinearConstraint> inequalities;
    
    //and keep track of AAind and rot...
    int[][] AAind, rot;
    int numRes;
    int numETerms;

    
    public PartialMinConstrainer(PairwiseEnergyMatrix eMatrix, CETMatrix cetm, boolean[][][][][][] sp,
            boolean[][][][][][][][][] tp, int AAind[][], int[][] rot) {
        this.eMatrix = eMatrix;
        this.cetm = cetm;
        
        splitFlags = sp;
        tripleFlags = tp;
        
        this.AAind = AAind;
        this.rot = rot;
        
        numRes = AAind.length;
        numETerms = numRes*(numRes+1)/2;//intra+shell, then pairwise
    }
    
    
    /*public PartialMinConstrainer(PairwiseEnergyMatrix eMatrix, PairwiseEnergyMatrix pairE) {
        this.eMatrix = eMatrix;
        this.pairE = pairE;
    }*/

    
    
    
    void addTupleEConstr(int... res){
        double coeffs[] = new double[numETerms];
        //intra+shell and pairwise energies limited to interactions in res
            
        for(int j : res)
            coeffs[j] = 1;
       
        //pairwise
        int count = numRes;
        for(int j=0; j<numRes; j++){
            for(int k=0; k<j; k++){
                if( RotamerSearch.arrayContains(res,j) && RotamerSearch.arrayContains(res,k) ){
                    coeffs[count] = 1;
                }
                count++;
            }
        }
        
        double minE = minTupleE(res);
        inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );
    }
    
    void addPairwiseEConstr(int res1, int res2){
        double coeffs[] = new double[numETerms];
       
        //pairwise
        int count = numRes;
        for(int j=0; j<numRes; j++){
            for(int k=0; k<j; k++){
                if( (res1==j&&res2==k) || (res1==k&&res2==j) ){
                    coeffs[count] = 1;
                }
                count++;
            }
        }
        
        double minE = getMinPairwiseE(res1,AAind[res1],rot[res1],res2,AAind[res2],rot[res2]);
        inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );
    }
    
    
    void addUndecidedTupleConstr(int... res){
        //Bound energy for the system consisting of all decided residues + 
        //and "undecided tuple" represented by res
        //this is designed to be easily evaluated for continuous sytems
        
        double coeffs[] = new double[numETerms];
        
        for(int j : res)
            coeffs[j] = 1;
       
        //pairwise
        int count = numRes;
        for(int j=0; j<numRes; j++){
            
            if(AAind[j].length==1)//j decided
                coeffs[j] = 1;
            
            for(int k=0; k<j; k++){
                //we count (j,k) iff each of j and k is either decided or in res
                if( (RotamerSearch.arrayContains(res,j)||AAind[j].length==1) && 
                        (RotamerSearch.arrayContains(res,k)||AAind[k].length==1) ){
                    coeffs[count] = 1;
                }
                count++;
            }
        }
        
        //DEBUG!!
        //don't add a constraint redundant with a previous one
        //we can check this in a more efficient way later if needed
        for(LinearConstraint lc : inequalities){
            boolean identical=true;
            for(int t=0; t<numETerms; t++){
                if(lc.getCoefficients().getEntry(t)!=coeffs[t]){
                    identical = false;
                    break;
                }
            }
            if(identical)
                return;//don't add constraint
        }
        //also don't add an all-zero constraint
        boolean allZero = true;
        for(int t=0; t<numETerms; t++){
            if(coeffs[t]!=0){
                allZero = false;
                break;
            }
        }
        if(allZero)
            return;
        
        
        
        double minE = minUndecidedTupleE(res);
        inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );
    }
    
    
    
    void addTupCombConstr(int maxRes){
        //add constraints on all combinations of energy terms with up to maxRes 
        //undecided residues involved (up to totMaxRes residues involved including decided)
        
        //in order to add each term just once, we'll create them in ascending order of term number
        //adding additional ones until no more are possible
        
        Deque<ETermComb> activeComb = new ArrayDeque<>();
        ArrayList<ETermComb> readyComb = new ArrayList<>();
        
        //we'll start with all fully decided energy terms
        double decidedCoeffs[] = new double[numETerms];
        int count = numRes;//counter for pairwise terms
        int totMaxRes = maxRes;
        HashSet<Integer> decidedRes = new HashSet<>();
        
        for(int j=0; j<numRes; j++){
            
            if(AAind[j].length==1){//j decided
                decidedCoeffs[j] = 1;
                totMaxRes++;
                decidedRes.add(j);
            }
            
            for(int k=0; k<j; k++){
                if(AAind[k].length==1&&AAind[j].length==1)
                    decidedCoeffs[count] = 1;
                count++;
            }
        }
        
        //put on stack
        activeComb.addFirst(new ETermComb(decidedCoeffs,decidedRes));
        
        
        while(!activeComb.isEmpty()){
            ETermComb curComb = activeComb.pollFirst();
            
            if(!curComb.resInvolved.isEmpty())//no point adding empty constrs
                readyComb.add(curComb);
            
            int firstTermOption = numETerms;//first of the E-terms we can assign
            //we can assign any term from the section of 0s at the end of curComb.coeffs
            //unless they add more residues and we're maxed out
            //(EDIT: we skip decided coefficients)
            while(curComb.coeffs[firstTermOption-1]==0||decidedCoeffs[firstTermOption-1]==1){
                firstTermOption--;
                if(firstTermOption==0)
                    break;
            }
            
            for(int newTerm=firstTermOption; newTerm<numETerms; newTerm++){
                
                if(decidedCoeffs[newTerm]==1)//SKIP
                    continue;
                
                HashSet<Integer> termResInvolved = resInvolvedInTerm(newTerm);
                                
                HashSet<Integer> newCombResInvolved = new HashSet<>();
                for(int ri : curComb.resInvolved)
                    newCombResInvolved.add(ri);
                for(int ri : termResInvolved)
                    newCombResInvolved.add(ri);
                
                if(newCombResInvolved.size() > totMaxRes)
                    //can't use this combination...too many residues involved
                    continue;
                
                double[] newCoeffs = curComb.coeffs.clone();
                newCoeffs[newTerm] = 1;
                
                
                activeComb.addLast(new ETermComb(newCoeffs,newCombResInvolved));
            }
        }
        
        for(ETermComb comb : readyComb){
            double minE = minTermCombE(comb);//make this similar to tuple E...
            inequalities.add( new LinearConstraint(comb.coeffs,Relationship.GEQ,minE) );
            
            //if there are some decided res we also need a version with the decided coeffs stripped out...
            if(!decidedRes.isEmpty()&&!(comb.resInvolved.size()==decidedRes.size())){
                double nodecCoeffs[] = comb.coeffs.clone();
                for(int c=0; c<numETerms; c++){
                    if(decidedCoeffs[c]>0)
                        nodecCoeffs[c] = 0;
                }
                ETermComb nodec = new ETermComb(nodecCoeffs);
                double nodecMinE = minTermCombE(nodec);
                inequalities.add( new LinearConstraint(nodec.coeffs,Relationship.GEQ,nodecMinE) );
            }
        }
    }
    
    
    HashSet<Integer> resInvolvedInTerm(int termNum){
        //for a given energy term, what are the residue involved?
        HashSet<Integer> ans = new HashSet<>();
        
        if(termNum<numRes)//intra+shell term
            ans.add(termNum);
        else {
            int count = numRes;
            for(int j=0; j<numRes; j++){
                if(termNum-count<j){
                    ans.add(j);
                    ans.add(termNum-count);
                    break;
                }
                count += j;
            }
        }
        
        return ans;
    }
    
    private class ETermComb {
        //combination of energy terms, along with the residues involved
        double[] coeffs;
        HashSet<Integer> resInvolved;

        public ETermComb(double[] coeffs, HashSet<Integer> resInvolved) {
            this.coeffs = coeffs;
            this.resInvolved = resInvolved;
        }
        
        public ETermComb(double[] coeffs) {
            //generate resInvolved
            this.coeffs = coeffs;
            resInvolved = new HashSet<>();
            for(int c=0; c<coeffs.length; c++){
                if(coeffs[c]!=0){
                    for(int ri : resInvolvedInTerm(c))
                        resInvolved.add(ri);
                }
            }
        }
    }
    
    
    double minTermCombE(ETermComb comb){
        
        ArrayList<RotTuple> possibleTups = new ArrayList<>();
        possibleTups.add(new RotTuple());//start with empty tuple
        
        //add all options at each residue in res
        for(int r : comb.resInvolved){
            ArrayList<RotTuple> newTups = new ArrayList<>();
            for(int j=0; j<AAind[r].length; j++){//assumed to be the same as rot.length
                for(RotTuple rt : possibleTups){
                    if(!rt.newRotIncompatible(r,AAind[r][j],rot[r][j])){
                        RotTuple bigTup = rt.copy();
                        bigTup.addRot(r,AAind[r][j],rot[r][j]);
                        newTups.add(bigTup);
                    }
                }
            }
            possibleTups = newTups;
        }
        
        double ans = Double.POSITIVE_INFINITY;
        for(RotTuple rt : possibleTups)
            ans = Math.min(ans,rt.energy(comb.coeffs));
        
        return ans;
    }
    
    
    
    
    void addDominantConstr(int... res){
        
        double coeffs[] = new double[numETerms];
            
        //coeffs applies to all intra+shell energies, plus all pairwise interactions involving res       
        
        //intra+shell
        for(int j=0; j<numRes; j++){
            if( scaleResDom && !RotamerSearch.arrayContains(res,j) )
                coeffs[j] = 1./(numRes-res.length);
            else
                coeffs[j] = 1;
        }

        //pairwise
        int count = numRes;
        for(int j=0; j<numRes; j++){
            for(int k=0; k<j; k++){
                if( RotamerSearch.arrayContains(res,j) || RotamerSearch.arrayContains(res,k) ){
                    coeffs[count] = 1;
                }
                count++;
            }
        }


        //add dominant constraint
        double minE = minResDominantE(res);
        inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );
    }
    
    //the one-res STAR constr
    //STAR STAR STAR
    void addOneResStarConstr(int res){
        double coeffsStar[] = new double[numETerms];

        //star applies to res intra+shell plus pairwise interactions of lower-numbered residues with it
        
        coeffsStar[res] = 1;

        int count = numRes;
        for(int j=0; j<numRes; j++){
            for(int k=0; k<j; k++){
                if(k==res||j==res){
                    if(k<res||j<res)//the other residue is <res
                        coeffsStar[count] = 1;
                }
                count++;
            }
        }
        
        double minStarE = minOneResStarE(res, AAind, rot);
        inequalities.add( new LinearConstraint(coeffsStar,Relationship.GEQ,minStarE) );
    }
    //we can try variations of this!  ordering up instead of down, by number of rots etc
    //current ordering has advantage for fixed ordering (group interactions with fixed res)
    
    //OK basically we want to lower-bound the energy for a conformational space
    //by enforcing lower-bound energies for parts of the system
    //like E_intra+shell (i) + intra+shell(j) + pairwise(i,j)
    //allowed rots at position res are (AAind[res][j],rot[res][j]) for each j
    double getBound(){
        //given conformation, compute bound
        
        //all these constraints are linear in the pairwise + intra energies
        //so we have an LP problem to minimize total energy
        //Let's use SimplexSolver
        
                        
        inequalities = new ArrayList<>();
        
        if(doTupComb)
            addTupCombConstr(2);
        else {
            //basic intra+shell, pairwise, and pair-energy constraints
            //should be enough to bound solution from below
            for(int res=0; res<numRes; res++){
                addTupleEConstr(res);
                for(int res2=0; res2<res; res2++){
                    addTupleEConstr(res,res2);
                    addPairwiseEConstr(res,res2);
                }
            }


            if(decidedBased){
                //good for cont flex because tupleE, undecidedTuple constr support that

                addUndecidedTupleConstr();//tuple constr for decided res

                for(int res=0; res<numRes; res++){
                    if(AAind[res].length>1){//undecided
                        addUndecidedTupleConstr(res);

                        if(do2Res){
                            for(int res2=0; res2<res; res2++){
                                if(AAind[res2].length>1){
                                    addUndecidedTupleConstr(res,res2);

                                    if(do3Res){
                                        for(int res3=0; res3<res2; res3++){
                                            if(AAind[res3].length>1){
                                                addUndecidedTupleConstr(res,res2,res3);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                //add star and 1-res dominant constraints
                //2 and 3 optionally too

                for(int res=0; res<numRes; res++){
                    addOneResStarConstr(res);
                    addDominantConstr(res);

                    if(do2Res){
                        for(int res2=0; res2<res; res2++){
                            addDominantConstr(res,res2);
                            if(do3Res){
                                for(int res3=0; res3<res2; res3++)
                                    addDominantConstr(res,res2,res3);
                            }
                        }
                    }
                }
            }
        }
        
        
        //any of these being infinite indicates clash
        for(LinearConstraint lc : inequalities){
            if(Double.isInfinite(lc.getValue()))
                return Double.POSITIVE_INFINITY;
        }
        
        
        //intra+shell
        /*for(int res=0; res<numRes; res++){
            double coeffs[] = new double[numETerms];
            coeffs[res] = 1;
            double minISE = getMinIntraAndShellE(res, AAind[res], rot[res]);
            
            inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minISE) );
        }
        
        int varCount = numRes;
        
        //pairwise and full pair E
        for(int res=0; res<numRes; res++){
            for(int res2=0; res2<res; res2++){
                
                //pairwise E
                double coeffs[] = new double[numETerms];
                coeffs[varCount] = 1;
                double minPairwiseE = getMinPairwiseE(res, AAind[res], rot[res], res2, AAind[res2], rot[res2]);
                
                if(Double.isInfinite(minPairwiseE))//if all pairs are pruned at some res pair, we can prune this node
                    return Double.POSITIVE_INFINITY;
                
                inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minPairwiseE) );
                
                //full-pair E
                double coeffsFP[] = new double[numETerms];
                coeffsFP[varCount] = 1;
                coeffsFP[res] = 1;
                coeffsFP[res2] = 1;
                
                //double minPairE = pairE.getPairwiseE(res, AAind[res], rot[res], res2, AAind[res2], rot[res2]);
                double minPairE = getPairMinE(res, AAind[res], rot[res], res2, AAind[res2], rot[res2]);
                
                inequalities.add( new LinearConstraint(coeffsFP,Relationship.GEQ,minPairE) );                
                varCount++;
            }
        }
        
        
        //ok now 1-res-dominant constr
        for(int res=0; res<numRes; res++){
            double coeffs[] = new double[numETerms];
            double coeffsStar[] = new double[numETerms];
            
            //coeffs applies to all intra+shell energies, plus pairwise interactions with res
            //star just to res intra+shell plus pairwise interactions of lower-numbered residues with it
            for(int j=0; j<numRes; j++){
                if(scaleResDom&&j!=res)
                    coeffs[j] = 1./(numRes-1);
                else
                    coeffs[j] = 1;
            }
            
            coeffsStar[res] = 1;
            
            int count = numRes;
            for(int j=0; j<numRes; j++){
                for(int k=0; k<j; k++){
                    if(k==res||j==res){
                        coeffs[count] = 1;
                        if(k<res||j<res)//the other residue is <res
                            coeffsStar[count] = 1;
                    }
                    count++;
                }
            }
            
            
            //add dominant constraint
            double minE = minOneResDominantE(res, AAind, rot);
            inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );
            
            
            //now add star constraint in the same way
            double minStarE = minOneResStarE(res, AAind, rot);
            inequalities.add( new LinearConstraint(coeffsStar,Relationship.GEQ,minStarE) );
        }
        
        
       
        //two-res-dominant and star
        if(do2Res){
            for(int res=0; res<numRes; res++){
                for(int res2=0; res2<res; res2++){
                    double coeffs[] = new double[numETerms];
                    double coeffsStar[] = new double[numETerms];

                    //coeffs applies to all intra+shell energies, plus pairwise interactions with res
                    //star just to res intra+shell plus pairwise interactions of lower-numbered residues with it
                    for(int j=0; j<numRes; j++){
                        if(scaleResDom&&j!=res&&j!=res2)
                            coeffs[j] = 1./(numRes-2);
                        else
                            coeffs[j] = 1;
                    }

                    coeffsStar[res] = 1;
                    coeffsStar[res2] = 1;

                    int count = numRes;
                    for(int j=0; j<numRes; j++){
                        for(int k=0; k<j; k++){
                            if(k==res||j==res||k==res2||j==res2){
                                coeffs[count] = 1;
                                if(k<res||j<res)//the other residue is <res
                                    coeffsStar[count] = 1;
                            }
                            count++;
                        }
                    }



                    //add dominant constraint
                    double minE = minTwoResDominantE(res, res2, AAind, rot);
                    inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );
 


                    
                }
            }
        }
        
        
        if(do3Res){
            for(int res=0; res<numRes; res++){
                for(int res2=0; res2<res; res2++){
                    for(int res3=0; res3<res2; res3++){
                        double coeffs[] = new double[numETerms];
                        double coeffsStar[] = new double[numETerms];

                        //coeffs applies to all intra+shell energies, plus pairwise interactions with res
                        //star just to res intra+shell plus pairwise interactions of lower-numbered residues with it
                        for(int j=0; j<numRes; j++){
                            if(scaleResDom&&j!=res&&j!=res2 && j!=res3)
                                coeffs[j] = 1./(numRes-3);
                            else
                                coeffs[j] = 1;
                        }

                        coeffsStar[res] = 1;
                        coeffsStar[res2] = 1;

                        int count = numRes;
                        for(int j=0; j<numRes; j++){
                            for(int k=0; k<j; k++){
                                if(k==res||j==res||k==res2||j==res2||j==res3||k==res3){
                                    coeffs[count] = 1;
                                    if(k<res||j<res)//the other residue is <res
                                        coeffsStar[count] = 1;
                                }
                                count++;
                            }
                        }



                        //add dominant constraint
                        double minE = minThreeResDominantE(res, res2, res3, AAind, rot);
                        inequalities.add( new LinearConstraint(coeffs,Relationship.GEQ,minE) );


                        //NOT USING STAR FOR NOW
                        //now add star constraint in the same way
                        
                    }
                }
            }
        }*/
        
        
        //could add higher tuples here...
        
        double fullECoeffs[] = new double[numETerms];
        Arrays.fill(fullECoeffs, 1);
        
        
        
        //DEBUG!!!!
        //making sure the GMEC is within the inequalities for the root node
        /*int GMECAA[] = ;
        int GMECRot[] = ;
        double GMECE[] = new double[numETerms];
        for(int res=0; res<numRes; res++)
            GMECE[res] = eMatrix.getIntraAndShellE(res, GMECAA[res], GMECRot[res]);
        
        int count = numRes;
            
        for(int res=0; res<numRes; res++){
            for(int res2=0; res2<res; res2++){
                GMECE[count] = eMatrix.getPairwiseE(res, GMECAA[res], GMECRot[res],
                        res2, GMECAA[res2], GMECRot[res2]);
                count++;
            }
        }*/
        
        //now go through each ineq and make sure GMEC satisfies them all (not strictly though)...
        
        //DEBUG!!!
        

        LinearObjectiveFunction fullE = new LinearObjectiveFunction(fullECoeffs,0);
        double bound = 0;
        try {
            SimplexSolver ss = new SimplexSolver(1e-3,10);
            ss.setMaxIterations(1000);
            bound = ss.optimize(fullE, inequalities, GoalType.MINIMIZE, false).getValue();
        }
        catch(UnboundedSolutionException e){
                        
            System.err.println(e.getMessage());

            
        }
        
        

        //double[] finalEs = opt.getOptimizationResponse().getSolution();
        //double bound = fullE.value(finalEs);
        return bound;
    }
   
    
    double minTupleE(int... res){
        //first version supports continuous flex through getRotTupleEnergy
        
        ArrayList<RotTuple> possibleTups = new ArrayList<>();
        possibleTups.add(new RotTuple());//start with empty tuple
        
        //add all options at each residue in res
        for(int r : res){
            ArrayList<RotTuple> newTups = new ArrayList<>();
            for(int j=0; j<AAind[r].length; j++){//assumed to be the same as rot.length
                for(RotTuple rt : possibleTups){
                    if(!rt.newRotIncompatible(r,AAind[r][j],rot[r][j])){
                        RotTuple bigTup = rt.copy();
                        bigTup.addRot(r,AAind[r][j],rot[r][j]);
                        newTups.add(bigTup);
                    }
                }
            }
            possibleTups = newTups;
        }
        
        double ans = Double.POSITIVE_INFINITY;
        for(RotTuple rt : possibleTups)
            ans = Math.min(ans,rt.energy());
        
        //DEBUG!!
        /*double check = 0;
        
        
        //supporting intra+shell or pair...kind of a cop-out but should do for now
        if(res.length==1)
            check = getMinIntraAndShellE(res[0], AAind[res[0]], rot[res[0]]);
        else if(res.length==2)
            check = getPairMinE(res[0], AAind[res[0]], rot[res[0]],
                    res[1], AAind[res[1]], rot[res[1]]);
        //else
        //    throw new RuntimeException("ERROR: MinTupleE doesn't support >=triples yet!");
        */
        
        return ans;
    }
    
    
    double minResDominantE(int... res){
        if(res.length==1)
            return minOneResDominantE(res[0]);
        else if(res.length==2)
            return minTwoResDominantE(res[0],res[1]);
        else if(res.length==3)
            return minThreeResDominantE(res[0],res[1],res[2]);
        else
            throw new RuntimeException("ERROR: res-dominant doesn't support >triples yet!");
    }
    
    
    double minUndecidedTupleE(int... res){
        //lower bound on undecided tuple (+ all decided res) energy
        //very similar to minTupleE but we start with the decided tuple and add onto that
        
        ArrayList<RotTuple> possibleTups = new ArrayList<>();
        //start with a tuple consisting of everything decided
        RotTuple decided = new RotTuple();
        for(int j=0; j<numRes; j++){
            if(AAind[j].length==1)//j is decided
                decided.addRot(j, AAind[j][0], rot[j][0]);
        }
        possibleTups.add(decided);
        
        
        //add all options at each residue in res
        for(int r : res){
            
            //make sure r isn't decided, that makes no sense
            if(AAind[r].length==1)
                throw new RuntimeException("ERROR: Residue "+r+" is decided but included in undecided tuple");
            
            ArrayList<RotTuple> newTups = new ArrayList<>();
            for(int j=0; j<AAind[r].length; j++){//assumed to be the same as rot.length
                for(RotTuple rt : possibleTups){
                    if(!rt.newRotIncompatible(r,AAind[r][j],rot[r][j])){
                        RotTuple bigTup = rt.copy();
                        bigTup.addRot(r,AAind[r][j],rot[r][j]);
                        newTups.add(bigTup);
                    }
                }
            }
            possibleTups = newTups;
        }
        
        double ans = Double.POSITIVE_INFINITY;
        for(RotTuple rt : possibleTups)
            ans = Math.min(ans,rt.energy());
        
        return ans;
    }
    
    //need func to do energy for rot tuple...
    private class RotTuple {
        //a tuple of rotamers at different residues
        //these are corresponding lists of res, AA, rot
        ArrayList<Integer> res = new ArrayList<>(), AA = new ArrayList<>(), rot = new ArrayList<>();
        
        void addRot(int p, int a, int r){
            res.add(p);
            AA.add(a);
            rot.add(r);
        }
        
        boolean newRotIncompatible(int p, int a, int r){
            //Are any combinations of (p,a,r) with part of this RotTuple pruned?
            //we assume (p,a,r) is not pruned by itself
            
            for(int j=0; j<res.size(); j++){
                //check pair pruning
                if(isPairPruned(p,a,r,res.get(j),AA.get(j),rot.get(j)))
                    return true;
                
                if(tripleFlags!=null){//also check triple pruning
                    for(int j2=0; j2<j; j2++){
                        if(isTriplePruned(p,a,r,res.get(j),AA.get(j),rot.get(j),
                                res.get(j2),AA.get(j2),rot.get(j2))){
                            return true;
                        }
                    }
                }
            }
            
            return false;//found no pruning, so compatible
        }
        
        RotTuple copy(){
            RotTuple ans = new RotTuple();
            for(int j=0; j<res.size(); j++)
                ans.addRot(res.get(j), AA.get(j), rot.get(j));
            return ans;
        }
        
        
        double energy(){
            //all intra+shell for the rotamers, and pairwise between them
            //not including shell-shell!
            
            double ans = 0;
            
            if(cetm==null){
                for(int j=0; j<res.size(); j++){
                    ans += eMatrix.getIntraAndShellE(res.get(j), AA.get(j), rot.get(j));
                    for(int k=0; k<j; k++){
                        ans += eMatrix.getPairwiseE( res.get(j), AA.get(j), rot.get(j),
                                res.get(k), AA.get(k), rot.get(k) );
                    }
                }
            }
            else {
                
                ContETerm terms[] = new ContETerm[res.size()*(res.size()+1)/2];
                
                int count = 0;
                
                for(int j=0; j<res.size(); j++){
                    terms[count] = cetm.intraAndShellBounds[res.get(j)][AA.get(j)][rot.get(j)];
                    count++;
                    
                    for(int k=0; k<j; k++){
                        terms[count] = cetm.pairwiseBounds[res.get(j)][AA.get(j)][rot.get(j)][res.get(k)][AA.get(k)][rot.get(k)];
                        count++;
                    }
                }
                
                
                CETObjFunction of = cetm.getObjFunc(terms, !contOnly, null);
                CCDMinimizer ccdMin = new CCDMinimizer(of,false);
            

                DoubleMatrix1D optDOFs = ccdMin.minimize();
                if(optDOFs==null)//for ellipse bounds, if the highest ellipses don't all intersect together,
                    //we can exclude the conformations involving them (set bound to inf)
                    return Double.POSITIVE_INFINITY;

                return of.getValue( optDOFs );


                /*if(es.useSVE){
                    //cof minimized m...revert to unminimized state
                    m.updateCoordinates();
                    m.revertPertParamsToCurState();
                }*/
            }
             
            return ans;
        }
        
        
        double energy(double[] coeffs){
            //energy where terms are waited based on coeffs (coefficients for all energy terms)
            
            HashMap<Integer,Integer> flex2Tup = new HashMap<>();
            //map overall flexible residue indices to indices in this tuples
            //(basically reverse lookup for res)
            for(int j=0; j<res.size(); j++)
                flex2Tup.put(res.get(j),j);
            
            double ans = 0;
            
            ContETerm terms[] = null;
            if(cetm!=null){
                int termCount = 0;
                for(int t=0; t<coeffs.length; t++){
                    if(coeffs[t]==1)
                        termCount++;
                    else if(coeffs[t]!=0)
                        throw new RuntimeException("ERROR: Currently just minimizing unweighted sums of EPIC terms!");
                }
                terms = new ContETerm[termCount];
            }
            
            int termCtr = 0;
            
            for(int j=0; j<numRes; j++){
                if(flex2Tup.containsKey(j)){
                    int ji = flex2Tup.get(j);
                    if(cetm==null)
                        ans += coeffs[j] * eMatrix.getIntraAndShellE(res.get(ji), AA.get(ji), rot.get(ji));
                    else if(coeffs[j]==1){
                        terms[termCtr] = cetm.intraAndShellBounds[res.get(ji)][AA.get(ji)][rot.get(ji)];
                        termCtr++;
                    }
                }
            }
            
            int count = numRes;
            for(int j=0; j<numRes; j++){
                for(int k=0; k<j; k++){
                    if(flex2Tup.containsKey(j)&flex2Tup.containsKey(k)){
                        int ji = flex2Tup.get(j);
                        int ki = flex2Tup.get(k);
                        if(cetm==null){
                            ans += coeffs[count] * eMatrix.getPairwiseE( res.get(ji), AA.get(ji), rot.get(ji),
                              res.get(ki), AA.get(ki), rot.get(ki) );
                        }
                        else if(coeffs[count]==1){
                            terms[termCtr] = cetm.pairwiseBounds[res.get(ji)][AA.get(ji)][rot.get(ji)][res.get(ki)][AA.get(ki)][rot.get(ki)];
                            termCtr++;
                        }
                    }
                    count++;
                }
            }
            
            if(cetm==null)
                return ans;
            else {
                CETObjFunction of = cetm.getObjFunc(terms, !contOnly, null);
                CCDMinimizer ccdMin = new CCDMinimizer(of,false);
            

                DoubleMatrix1D optDOFs = ccdMin.minimize();
                if(optDOFs==null)//for ellipse bounds, if the highest ellipses don't all intersect together,
                    //we can exclude the conformations involving them (set bound to inf)
                    return Double.POSITIVE_INFINITY;

                return of.getValue( optDOFs );
            }
        }

    }
    
    

    double getMinIntraAndShellE(/*eMatrix.getIntraAndShellE(*/int res, int AAind[], int rot[]){
        //get minimum intra and shell energy for the given residue, with specified AA and rot options
        double ans = Double.POSITIVE_INFINITY;
 
        for(int j=0; j<AAind.length; j++)//assumed to be the same as rot.length
            ans = Math.min(ans,eMatrix.getIntraAndShellE(res, AAind[j], rot[j]));
        //and no continuous contribution is expected
        //WOULD BE GOOD TO SCREEN FOR PRUNED ROTS TOO
        
        return ans;
    }
    
    double getMinPairwiseE(int res1, int[] AAind1, int[] rot1, int res2, int[] AAind2, int[] rot2){
        double ans = Double.POSITIVE_INFINITY;
        
        for(int j=0; j<AAind1.length; j++){//assumed to be the same as rot1.length
            for(int k=0; k<AAind2.length; k++){//same as rot2.length
                if(!isPairPruned(res1, AAind1[j], rot1[j], res2, AAind2[k], rot2[k]))
                    ans = Math.min(ans,eMatrix.getPairwiseE(res1, AAind1[j], rot1[j], res2, AAind2[k], rot2[k]));
            }
        }
        //and no continuous contribution is expected
        
        return ans;
    }
    
    
    //pruning is even more key here...
    //another easy linear constraint is minimized energy of (unassigned res) interacting with all assigned res
    double getPairMinE(int res1, int[] AAind1, int[] rot1, int res2, int[] AAind2, int[] rot2){
        double ans = Double.POSITIVE_INFINITY;
        
        for(int j=0; j<AAind1.length; j++){//assumed to be the same as rot1.length
            for(int k=0; k<AAind2.length; k++){//same as rot2.length
                
                if(!isPairPruned(res1, AAind1[j], rot1[j], res2, AAind2[k], rot2[k])){
                    double pairE = eMatrix.getIntraAndShellE(res1, AAind1[j], rot1[j])
                        + eMatrix.getIntraAndShellE(res2, AAind2[k], rot2[k])
                        + eMatrix.getPairwiseE(res1, AAind1[j], rot1[j], res2, AAind2[k], rot2[k]);

                    if(cetm!=null)//possibility of wanting to cache this...
                        throw new RuntimeException("CETM NOT SUPPORTED YET");

                    ans = Math.min(ans,pairE);
                }
            }
        }
        
        return ans;
    }
    
    
    double minOneResDominantE(int res){
        //lower bound on the energy of res plus its interactions with other residues (including their intra+shell E's)
        double ans = Double.POSITIVE_INFINITY;
        
        if(cetm!=null)//possibility of wanting to cache this...
            throw new RuntimeException("CETM NOT SUPPORTED YET");
        
        int numRes = AAind.length;
        
        //loop through rotamers of res
        for(int j=0; j<AAind[res].length; j++){
            
            double rotE = eMatrix.getIntraAndShellE(res, AAind[res][j], rot[res][j]);
            
            for(int res2=0; res2<numRes; res2++){//go through all other residues
                if(res2!=res){
                    
                    double minInteraction = Double.POSITIVE_INFINITY;
                    
                    for(int k=0; k<AAind[res2].length; k++){
                        
                        if(!isPairPruned(res, AAind[res][j], rot[res][j], res2, AAind[res2][k], rot[res2][k])){
                            if(scaleResDom){
                                minInteraction = Math.min( minInteraction,
                                    eMatrix.getIntraAndShellE(res2, AAind[res2][k], rot[res2][k]) / (numRes-1) + 
                                    eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res2, AAind[res2][k], rot[res2][k]) );
                            }
                            else{
                                minInteraction = Math.min( minInteraction,
                                    eMatrix.getIntraAndShellE(res2, AAind[res2][k], rot[res2][k]) + 
                                    eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res2, AAind[res2][k], rot[res2][k]) );
                            }
                        }
                    }
                    
                    rotE += minInteraction;
                }
            }
            
            ans = Math.min(ans,rotE);
        }
        
        return ans;
    }
    
    
    double minOneResStarE(int res, int[][] AAind, int[][] rot){
        //lower bound on the energy of res plus its interactions with lower-numbered residues
        //these are the lower bounds used in quick F-score
        //and the star E's for all residues add up to the full energy
        //this should ensure that the LP bound is no looser than quick F-score
        double ans = Double.POSITIVE_INFINITY;
        
        if(cetm!=null)//possibility of wanting to cache this...
            throw new RuntimeException("CETM NOT SUPPORTED YET");
        
        //loop through rotamers of res
        for(int j=0; j<AAind[res].length; j++){
            
            double rotE = eMatrix.getIntraAndShellE(res, AAind[res][j], rot[res][j]);
            
            for(int res2=0; res2<res; res2++){//go through all previous residues
                    
                double minInteraction = Double.POSITIVE_INFINITY;

                for(int k=0; k<AAind[res2].length; k++){
                    if(!isPairPruned(res, AAind[res][j], rot[res][j], res2, AAind[res2][k], rot[res2][k])){
                        minInteraction = Math.min( minInteraction,
                            eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res2, AAind[res2][k], rot[res2][k]) );
                    }
                }

                rotE += minInteraction;
            }
            
            ans = Math.min(ans,rotE);
        }
        
        return ans;
    }
    
    
    double minTwoResDominantE(int res, int res2){
        //lower bound on the energy of (res,res2) pair plus its interactions with other residues (including their intra+shell E's)
        double ans = Double.POSITIVE_INFINITY;
        
        if(cetm!=null)//possibility of wanting to cache this...
            throw new RuntimeException("CETM NOT SUPPORTED YET");
        
        int numRes = AAind.length;
        
        //loop through rotamer pairs for (res,res2)
        for(int j=0; j<AAind[res].length; j++){
            for(int j2=0; j2<AAind[res2].length; j2++){
                
                if(!isPairPruned(res, AAind[res][j], rot[res][j], res2, AAind[res2][j2], rot[res2][j2])){
            
                    double rotE = eMatrix.getIntraAndShellE(res, AAind[res][j], rot[res][j]) +
                            eMatrix.getIntraAndShellE(res2, AAind[res2][j2], rot[res2][j2]) +
                            eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res2, AAind[res2][j2], rot[res2][j2]);

                    for(int res3=0; res3<numRes; res3++){//go through all other residues
                        if(res3!=res && res3!=res2){

                            double minInteraction = Double.POSITIVE_INFINITY;

                            for(int k=0; k<AAind[res3].length; k++){
                                
                                if(!isTriplePruned(res, AAind[res][j], rot[res][j], res2, AAind[res2][j2], rot[res2][j2], res3, AAind[res3][k], rot[res3][k])){
                                    if(scaleResDom){
                                        minInteraction = Math.min( minInteraction,
                                            eMatrix.getIntraAndShellE(res3, AAind[res3][k], rot[res3][k]) / (numRes-2) + 
                                            eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res3, AAind[res3][k], rot[res3][k]) + 
                                            eMatrix.getPairwiseE(res2, AAind[res2][j2], rot[res2][j2], res3, AAind[res3][k], rot[res3][k]) );
                                    }
                                    else{
                                        minInteraction = Math.min( minInteraction,
                                            eMatrix.getIntraAndShellE(res3, AAind[res3][k], rot[res3][k]) + 
                                            eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res3, AAind[res3][k], rot[res3][k]) + 
                                            eMatrix.getPairwiseE(res2, AAind[res2][j2], rot[res2][j2], res3, AAind[res3][k], rot[res3][k]) );
                                    }
                                }
                            }

                            rotE += minInteraction;
                        }
                    }

                    ans = Math.min(ans,rotE);
                }
            }
        }
        
        return ans;
    }
    
    
    double minThreeResDominantE(int res, int res2, int res3){
        //lower bound on the energy of (res,res2,res3) triple plus its interactions with other residues (including their intra+shell E's)
        double ans = Double.POSITIVE_INFINITY;
        
        if(cetm!=null)//possibility of wanting to cache this...
            throw new RuntimeException("CETM NOT SUPPORTED YET");
        
        int numRes = AAind.length;
        
        //loop through rotamer pairs for (res,res2)
        for(int j=0; j<AAind[res].length; j++){
            for(int j2=0; j2<AAind[res2].length; j2++){
                for(int j3=0; j3<AAind[res3].length; j3++){
                    
                    if(!isTriplePruned( res, AAind[res][j], rot[res][j],
                            res2, AAind[res2][j2], rot[res2][j2], res3, AAind[res3][j3], rot[res3][j3] )){
            
                        double rotE = eMatrix.getIntraAndShellE(res, AAind[res][j], rot[res][j]) +
                                eMatrix.getIntraAndShellE(res2, AAind[res2][j2], rot[res2][j2]) +
                                eMatrix.getIntraAndShellE(res3, AAind[res3][j3], rot[res3][j3]) +
                                eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res2, AAind[res2][j2], rot[res2][j2]) +
                                eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res3, AAind[res3][j3], rot[res3][j3]) +
                                eMatrix.getPairwiseE(res2, AAind[res2][j2], rot[res2][j2], res3, AAind[res3][j3], rot[res3][j3]);

                        for(int res4=0; res4<numRes; res4++){//go through all other residues
                            if(res4!=res && res4!=res2 && res4!=res3){

                                double minInteraction = Double.POSITIVE_INFINITY;

                                for(int k=0; k<AAind[res4].length; k++){
                                    
                                    if(!isQuadPruned(res, AAind[res][j], rot[res][j],
                                                res2, AAind[res2][j2], rot[res2][j2], 
                                                res3, AAind[res3][j3], rot[res3][j3],
                                                res4, AAind[res4][k], rot[res4][k])){
                                        
                                        if(scaleResDom){
                                            minInteraction = Math.min( minInteraction,
                                                eMatrix.getIntraAndShellE(res4, AAind[res4][k], rot[res4][k]) / (numRes-3) + 
                                                eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res4, AAind[res4][k], rot[res4][k]) + 
                                                eMatrix.getPairwiseE(res2, AAind[res2][j2], rot[res2][j2], res4, AAind[res4][k], rot[res4][k]) +
                                                eMatrix.getPairwiseE(res3, AAind[res3][j3], rot[res3][j3], res4, AAind[res4][k], rot[res4][k]) );
                                        }
                                        else{
                                            minInteraction = Math.min( minInteraction,
                                                eMatrix.getIntraAndShellE(res4, AAind[res4][k], rot[res4][k]) + 
                                                eMatrix.getPairwiseE(res, AAind[res][j], rot[res][j], res4, AAind[res4][k], rot[res4][k]) + 
                                                eMatrix.getPairwiseE(res2, AAind[res2][j2], rot[res2][j2], res4, AAind[res4][k], rot[res4][k]) +
                                                eMatrix.getPairwiseE(res3, AAind[res3][j3], rot[res3][j3], res4, AAind[res4][k], rot[res4][k]) );
                                        }
                                    }
                                }

                                rotE += minInteraction;
                            }
                        }

                        ans = Math.min(ans,rotE);
                    }
                }
            }
        }
        
        return ans;
    }
    
    
    
    //checking tuple pruning
    //we use whatever pruned pair or triple info we have (may have null)
    boolean isPairPruned(int p1, int a1, int r1, int p2, int a2, int r2){
        if(splitFlags!=null)
            return splitFlags[p1][a1][r1][p2][a2][r2];

        return false;
    }
    
    boolean isTriplePruned(int p1, int a1, int r1, int p2, int a2, int r2, int p3, int a3, int r3){
        //if there is a pruned pair in this triple, then the triple is impossible
        if(isPairPruned(p1,a1,r1,p2,a2,r2))
            return true;
        if(isPairPruned(p1,a1,r1,p3,a3,r3))
            return true;
        if(isPairPruned(p2,a2,r2,p3,a3,r3))
            return true;
        
        
        
        if(tripleFlags!=null)
            return isPrunedTriple(p1,a1,r1,p2,a2,r2,p3,a3,r3);

        return false;
    }
    
    
    //copied from DEE.java
    boolean isPrunedTriple(int curPos1, int curAA1, int curRot1,
                int curPos2, int curAA2, int curRot2,
                int curPos3, int curAA3, int curRot3 ){
            //Checks if a given triple is pruned in tripleFlags
            //curPos1, curPos2, and curPos3 (positions among mutable residues of
            //residues 1, 2, and 3) are assumed to be all different (otherwise it's not a prunable triple)
            //If useTriples==true, calling this with pruned rotamers might yield a null-pointer exception
            //(because of the space-saving structure of tripleFlags)

            if( curPos1 > curPos2 ){
                if( curPos2 > curPos3 )
                    return tripleFlags[curPos1][curAA1][curRot1][curPos2][curAA2][curRot2][curPos3][curAA3][curRot3];
                else if ( curPos1 > curPos3 )
                    return tripleFlags[curPos1][curAA1][curRot1][curPos3][curAA3][curRot3][curPos2][curAA2][curRot2];
                else//index3 > index1 > index2
                    return tripleFlags[curPos3][curAA3][curRot3][curPos1][curAA1][curRot1][curPos2][curAA2][curRot2];
            }
            else{
                if( curPos1 > curPos3 )
                    return tripleFlags[curPos2][curAA2][curRot2][curPos1][curAA1][curRot1][curPos3][curAA3][curRot3];
                else if ( curPos2 > curPos3 )
                    return tripleFlags[curPos2][curAA2][curRot2][curPos3][curAA3][curRot3][curPos1][curAA1][curRot1];
                else
                    return tripleFlags[curPos3][curAA3][curRot3][curPos2][curAA2][curRot2][curPos1][curAA1][curRot1];
            }
    }
    
    
    boolean isQuadPruned(int p1, int a1, int r1, int p2, int a2, int r2, int p3, int a3, int r3,
            int p4, int a4, int r4){
        //currently don't have quad flags so just check if lower tuples are pruneds
        if(isTriplePruned(p1,a1,r1,p2,a2,r2,p3,a3,r3))
            return true;
        if(isTriplePruned(p1,a1,r1,p2,a2,r2,p4,a4,r4))
            return true;
        if(isTriplePruned(p1,a1,r1,p3,a3,r3,p4,a4,r4))
            return true;
        if(isTriplePruned(p2,a2,r2,p3,a3,r3,p4,a4,r4))
            return true;
        
        
        return false;
    }
    
    
    
    
    
    //This non-CETM version had a lot of error (probably minimization error)
    //also only supported one assignment per position and no pruning...
    /*double getBound(int[] AAind, int[] rot){
        //given conformation, compute bound
        
        //all these constraints are linear in the pairwise + intra energies
        //so we have an LP problem to minimize total energy
        //Let's use JOptimizer
        
        int numRes = AAind.length;
        int numETerms = numRes*(numRes+1)/2;//intra+shell, then pairwise
        
        OptimizationRequest or = new OptimizationRequest();
        
        
        //2 constraints for each pair, + 1 for each intra+shell energy
        LinearMultivariateRealFunction[] inequalities 
                        = new LinearMultivariateRealFunction[numRes*numRes];
        
        double initPt[] = new double[numETerms];//we'll just set this 1 above each constraint
        
        //intra+shell
        for(int res=0; res<numRes; res++){
            double coeffs[] = new double[numETerms];
            coeffs[res] = -1;
            double minISE = eMatrix.getIntraAndShellE(res, AAind[res], rot[res]);
            inequalities[res] = new LinearMultivariateRealFunction(coeffs,minISE);
            
            initPt[res] = minISE + 1;
        }
        
        int constrCount = numRes;
        int varCount = numRes;
        
        //pairwise and full pair E
        for(int res=0; res<numRes; res++){
            for(int res2=0; res2<res; res2++){
                
                //pairwise E
                double coeffs[] = new double[numETerms];
                coeffs[varCount] = -1;
                double minPairwiseE = eMatrix.getPairwiseE(res, AAind[res], rot[res], res2, AAind[res2], rot[res2]);
                inequalities[constrCount] = new LinearMultivariateRealFunction(coeffs,minPairwiseE);
                constrCount++;
                
                //full-pair E
                double coeffsFP[] = new double[numETerms];
                coeffsFP[varCount] = -1;
                coeffsFP[res] = -1;
                coeffsFP[res2] = -1;
                double minPairE = pairE.getPairwiseE(res, AAind[res], rot[res], res2, AAind[res2], rot[res2]);
                inequalities[constrCount] = new LinearMultivariateRealFunction(coeffsFP,minPairE);
                
                
                //initialization
                initPt[varCount] = minPairwiseE + 1;
                double minusInitPairE = inequalities[constrCount].value(initPt);
                if(minusInitPairE >= 0)
                    initPt[varCount] += minusInitPairE+1;
                
                
                constrCount++;
                varCount++;
            }
        }
        
        double fullECoeffs[] = new double[numETerms];
        Arrays.fill(fullECoeffs, 1);
        
        
        or.setFi(inequalities);
        or.setInitialPoint(initPt);
        
        or.setToleranceFeas(1.E-12);//as online
        or.setTolerance(1.E-12);
        //trying a little looser tolerance...
        //or.setToleranceFeas(0.001);
        //or.setTolerance(0.001);
        
        LinearMultivariateRealFunction fullE = new LinearMultivariateRealFunction(fullECoeffs,0);
        or.setF0(fullE);

        //optimization
        JOptimizer opt = new JOptimizer();
        opt.setOptimizationRequest(or);
        try{
            int returnCode = opt.optimize();
            System.out.println("JOptimizer return code: "+returnCode);
            //kinda too much output...
        }
        catch(Exception e){
            System.err.print("JOptimizer error: "+e.getMessage());
            e.printStackTrace();
            System.err.println("Outputting MaxResidMinimizerJopt to MRMJ.dat for inspection");
            KSParser.outputObject(this, "MRMJ.dat");
            System.exit(1);
        }

        double[] finalEs = opt.getOptimizationResponse().getSolution();
        double bound = fullE.value(finalEs);
        return bound;
    }*/
    
    
}
