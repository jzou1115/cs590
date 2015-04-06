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
//	CCDMinimizer.java
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
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.BitSet;


//This is a modular CCD minimizer
////with the objective function handling the specifics about the degrees of freedom
//(dihedrals, perturbation parameters, or whatever)

//It is meant to be able to reproduce the steps of the SimpleMinimizer
//but also to adjust the step lengths to allow faster
//A Quasi-Newton minimizer was also tried but didn't work very well (see revision 6 or earlier of ASDF svn)

public class CCDMinimizer {

    ObjectiveFunction objFcn;

    DoubleMatrix1D x;//Current values of degrees of freedom

    int numDOFs;

    static int numIter=30;
    static double EConvTol = 0.001;//convergence threshold
    static double numTol = 1e-6;
    static double GCTol = 1e-10;//tolerance for violation of non-box constraints

    double[] firstStep, lastStep;//Keep track of first and most recent (last) steps in each DOF

    boolean useCorners;//Indicates if we should search the corners of the voxel for
    //the best initial value or just use the middle

    //Constraint information
    DoubleMatrix1D DOFmin, DOFmax;

    DoubleMatrix1D singleInitVal = null;//If this isn't null,
    //then instead of normal initial value calculation we just try to minimize from this initial value
    //we return null from minimize() if it's out of range

    
    
    //The rest of the fields are pretty special-purpose (not needed for just minimizing energy over a box, etc.)
    
    //Support for non-box constraints, implemented using the GenCoord class
    //These non-box constraints need to be set before performing minimization if they are needed
    GenCoord nonBoxConstrGC[] = new GenCoord[0];
    int nonBoxConstrIndices[][];//Indices in x of DOFs that the GenCoords in nonBoxConstrGC operate on
    double GCmin[];
    double GCmax[];//The GenCoord nonBoxConstrGC[a] will be constrained to the range [ GCmin[a], GCmax[a] ]
    int nonBoxConstrAffectingDOF[][];//nonBoxConstrAffectingDOF[a] lists the non-box constraints affecting DOF a of x
    //(listed as indices in nonBoxConstrGC)

    GenCoord rescalingGC = null;//If this isn't null, at the end of each iteration, we will
    //scale x by rescalingGC(x)
    //used for minimization over level sets to keep the scaling factor from being too large or small
    int rescalingIndices[] = null;//Indices for rescalingGC to operate on
    //(can be kept null if rescalingGC doesn't use indices)

    boolean banZero = false;//do not allow an all-zero DOF vector (used for minimization over level sets)

    boolean useRandMinCheck = false;//do a randomized check of the minimum to make sure there aren't
    //in-range, lower-energy points nearby


    boolean jumpOOR = false;//jump out-of-range values due to GC's (intended for level-set minimization
    //or more generally for minimizing a tractable function over a disconnected subset of its domain)
    double jumpOORDelta = 0.05;//jumpOOR parameter defaults are for angles in radians
    double jumpOORNumSteps = 126;//enough to go all the way around


    double minTime;//time for most recent minimization (ms)

    public CCDMinimizer( ObjectiveFunction ofn, boolean useCorners ){

        objFcn = ofn;
        numDOFs = ofn.getNumDOFs();
        nonBoxConstrAffectingDOF = new int[numDOFs][0];
        //No constraints specified until minimize() is called
    }

    public DoubleMatrix1D minimize() {

        long minStartTime = System.currentTimeMillis();
        
        //First figure out the constraints and initial values
        //(since the minimizer might be used for many rotameric states, etc.,
        //this can't be done in the constructor)
        DoubleMatrix1D constr[] = objFcn.getConstraints();
        DOFmin = constr[0];
        DOFmax = constr[1];
        
        if(!compInitVals())//Initialize x
            return null;//No initial values found (this is currently only for IBEX)

        firstStep = new double[numDOFs];
        lastStep = new double[numDOFs];
        Arrays.fill(firstStep, 1);
        Arrays.fill(lastStep, 1);


        double oldE = Double.POSITIVE_INFINITY;

        recordCoordForFBFB();                               
        
        objFcn.setDOFs(x);
        
        
        //DEBUG!!!
        //new KSParser().saveMolecule( ((ContSCObjFunction)objFcn).m, "PREMIN2.pdb", ((ContSCObjFunction)objFcn).efunc.getEnergy() );

        
        for(int iter=0; iter<numIter; iter++) {

            double E = objFcn.getValue(x);
            
            if(Double.isInfinite(E)&&Double.isInfinite(oldE)){
                
                if(useRandMinCheck)
                    E = doRandMinCheck(E);

                if(Double.isInfinite(E))//if useRandMinCheck, make sure the check was passed
                    break;
            }
            //Infinite energy two steps in a row implies we're not getting out
            //Could also mean we started at infinity, but due to the initial value computation setup,
            //this should only happen for a uniformly infinite objective function

            //System.out.println("Iteration: "+iter+" Energy: "+E);

            //Considered to be converged if energy improvement is small enough
            if( oldE - E < EConvTol ){//or if numerical error is causing the energy to rise
                //(for EConvTol small enough, or for very rugged regions like clashes, the rise may exceed EConvTol)

                if(useRandMinCheck)
                    E = doRandMinCheck(E);

                if( oldE - E < EConvTol )//if using RandMinCheck, make sure it is passed
                    break;
            }


            oldE = E;

            for(int dof=0; dof<numDOFs; dof++){

                //Don't minimize for this DOF if it is constrained to a single value
                if( DOFmax.get(dof) == DOFmin.get(dof) )
                    continue;

                //Do line search for this dof

                double dof_base = x.get(dof);
                double curVal = objFcn.getValForDOF(dof,dof_base);//Get value for DOF, shifted up or down as needed
                double step = getStepSize(dof,iter);//Step size for the given degree of freedom, adjusted for iteration


                //Values up or down a step for this DOF: initialized to inf (indicating infeasibility)
                //in case the DOF value is out of range
                double upVal = Double.POSITIVE_INFINITY;
                double downVal = Double.POSITIVE_INFINITY;

                if( ! isOutOfRange(dof_base+step,dof) )
                    upVal = objFcn.getValForDOF(dof,dof_base+step);
                if( ! isOutOfRange(dof_base-step,dof) )
                    downVal = objFcn.getValForDOF(dof,dof_base-step);

                //quadratic approximation of obj function: curVal + a*(dof_val-dof_base)^2 + b*(dof_val-dof_base);
                //a*step^2 + b*step = upVal - curVal;
                //a*step^2 - b*step = downVal-curVal;
                //


                double a = (upVal+downVal-2*curVal)/(2*step*step);
                double estmin=0;

                
                //TRYING DIFFERENT VERSIONS
                if( ( a <= 0 ) || Double.isNaN(a) || Double.isInfinite(a) ){
                    //negative a is not a good sign...use the lesser of upVal, downVal
                    //infinite or nan a means we're hitting a constraint or impossible conformation
                    if( upVal < downVal )
                        estmin = dof_base+step;
                    else
                        estmin = dof_base-step;
                }
                
                
               /* if(Double.isInfinite(upVal)){//near upper constraint...need to be able to go right up to edge
                    //or away from edge as appropriate
                    if(downVal<curVal)
                        estmin = dof_base-step;
                    else
                        estmin = dof_base+step;
                }
                else if(Double.isInfinite(downVal)){//near lower constraint
                    if(upVal<curVal)
                        estmin = dof_base+step;
                    else
                        estmin = dof_base-step;
                }
                else if( ( a <= 0 ) || Double.isNaN(a) ){
                    //negative a is not a good sign...use the lesser of upVal, downVal
                    //infinite or nan a means we're hitting a constraint or impossible conformation
                    if( upVal < downVal )
                        estmin = dof_base+step;
                    else
                        estmin = dof_base-step;
                } */
                else{
                    double b = (upVal-curVal)/step-a*step;
                    estmin = dof_base-b/(2*a);//2*a*(dof_val-dof_base)+b=0 here, with positive a
                }


                if( isOutOfRange(estmin,dof) )
                    estmin = getEdgeDOFVal(estmin,dof);


                //if( Math.abs(estmin-dof_base) < numTol )//This can happen if both upVal and downVal are infinite (perhaps due to a loop closure failure)
                //    break;

                double estminVal = objFcn.getValForDOF(dof,estmin);
                double estminValOld = curVal;

                if(estminVal < curVal){
                    while(true) {//Can break either on hitting a constraint or on estminVal starting to increase
                        estmin = dof_base + 2*(estmin-dof_base);

                        if( isOutOfRange(estmin,dof) ){
                            double edge = getEdgeDOFVal(estmin,dof);

                            double edgeVal = objFcn.getValForDOF(dof,edge);

                            if(edgeVal<estminVal)
                                x.set(dof, edge);
                            else
                                x.set(dof, dof_base+0.5*(estmin-dof_base) );

                            break;
                        }

                        estminValOld = estminVal;
                        estminVal = objFcn.getValForDOF(dof,estmin);

                        if( !(estminVal < estminValOld) ){//No improvement in the last step
                            x.set(dof,dof_base+0.5*(estmin-dof_base));
                            break;
                        }
                    }
                }
                else if(estminVal>curVal+numTol) {//need to backtrack
                    //won't hit a constraint 
                     while(true) {//Can break on estminVal starting to increase, or decreasing negligibly
                        estmin = dof_base + 0.5*(estmin-dof_base);

                        estminValOld = estminVal;
                        estminVal = objFcn.getValForDOF(dof,estmin);

                        if( estminValOld < estminVal + numTol ){//No significant improvement in the last step
                            if(estminValOld<curVal)//have improvement over curVal at least
                                x.set(dof,dof_base+2*(estmin-dof_base));
                            break;
                        }
                    }
                }


                lastStep[dof] = x.get(dof) - dof_base;
                if(iter==0)
                    firstStep[dof] = lastStep[dof];


                //RIPPLE JUMPER!!
                //search up and down some, with the goal of being robust for rugged surfaces
                //optimized for dihedrals...
                double downRipple = x.get(dof)-1;
                double upRipple = x.get(dof)+1;
                curVal = objFcn.getValForDOF(dof, x.get(dof));
                if( ! isOutOfRange( downRipple, dof) ){
                    if( objFcn.getValForDOF(dof,downRipple) < curVal )
                        x.set(dof, downRipple);
                }
                if( ! isOutOfRange( upRipple, dof ) ){
                    if( objFcn.getValForDOF(dof,upRipple) < curVal )
                        x.set(dof, upRipple);
                }
                //END RIPPLE JUMPER

                objFcn.setDOF(dof, x.get(dof));
                
                recordCoordForFBFB();
            }
            
            if(rescalingGC != null)
                rescaleValues();
        }

        minTime = System.currentTimeMillis() - minStartTime;
        return x;

    }

    


    public double getStepSize(int dof, int iter){//Get a good step size for a DOF
        //for initial value checking in CCD, in iteration number iter

        double initStepSize = objFcn.getInitStepSize(dof);
        
        //trying to make it adaptive (based on historical steps if possible; else on step #)

        if( Math.abs(lastStep[dof]) > numTol && Math.abs(firstStep[dof]) > numTol )
            return initStepSize * Math.abs(lastStep[dof]/firstStep[dof]);
        else
            return initStepSize / Math.pow(iter + 1, 3);
    }

    public double getDOFTol(int dof, int iter){//Get a good error tolerance (absolute error in the DOF value) for
        //a given iteration of CCD
        return getStepSize(dof,iter)*4;//This seems to work pretty well
    }


    protected boolean isOutOfRange(double val, int dof) {//Is the value out of range for the degree of freedom?

        if( val<DOFmin.get(dof) || val>DOFmax.get(dof) )//regular box constraints
            return true;
        else if( nonBoxConstrAffectingDOF[dof].length > 0 ){//handling non-box constraints

            double curValForDOF = x.get(dof);
            x.set(dof, val);//temporarily set x for constraint checking

            if(banZero){
                if(x.zDotProduct(x)==0){
                    x.set(dof, curValForDOF);
                    return true;
                }
            }

            for( int nbc : nonBoxConstrAffectingDOF[dof] ) {
                double GCVal = nonBoxConstrGC[nbc].eval(x, nonBoxConstrIndices[nbc]);
                if( GCVal > GCmax[nbc] + GCTol || GCVal < GCmin[nbc] - GCTol ){//Out of range for non-box constraint
                    //We let values within a small tolerance of the border slide, because when we call getEdgeDOFVal we may
                    //be trying to get something right on the border that will be slightly off when the coordinate is evaluated
                    x.set(dof, curValForDOF);
                    return true;
                }
            }

            x.set(dof, curValForDOF);
            //set x back because we may not use this value
            //as an accepted minimization step
        }

        if(banZero){
            if(x.zDotProduct(x)==0)
                return true;
        }

        return false;
    }

    double getEdgeDOFVal(double oorVal, int dof){//checkx true by default
        //(but should be set to false if we're trying to get x in range using getEdgeDOFVal)
        return getEdgeDOFVal(oorVal, dof, true, 0);
    }

    double getEdgeDOFVal(double oorVal, int dof, boolean checkx){
        //we start at recursion 0
        return getEdgeDOFVal(oorVal, dof, checkx, 0);
    }



    double getEdgeDOFVal(double oorVal, int dof, boolean checkx, int numRecursions) {
        //Given the out-of-range value oorVal for DOF dof, return the
        //nearest value of dof that is at the edge of the feasible region
        //We assume the other DOFs are set correctly in x
        //and that some value of dof will be in range given the other DOFs
        double edge = oorVal;

        if(checkx){//check if x is in range...during minimization it should be (considering substituting oorVal into x)
            GCTol = 0.05;//don't be too strict on GCs
            if(isOutOfRange(x.get(dof),dof))//not a problem if trying to get initial values...
                System.out.println("x="+x.get(dof)+" out of range for getEdgeDOFVal.  DOF min: "
                        +DOFmin.get(dof)+" max: "+DOFmax.get(dof));
            GCTol = 1e-10;
        }

        if(banZero)
            edge = enforceZeroBan(edge,dof);


        //GenCoord handling.  
        //We will assume that we can move the point, fixing the other DOFs, to satisfy one constraint at a time
        //This should work e.g. if the nearest feasible point is in the same direction for all constraints (as is expected,
        //and is guaranteed if the intersection of the feasible region with the line we're searching on is connected)
        if( nonBoxConstrAffectingDOF[dof].length > 0 ){

            for( int nbc : nonBoxConstrAffectingDOF[dof] ){

                if(banZero)
                    edge = enforceZeroBan(edge,dof);

                edge = nonBoxConstrGC[nbc].getNearestInRangeDOFVal( edge,
                        GCmin[nbc], GCmax[nbc], x, dof, nonBoxConstrIndices[nbc] );
            }
        }

        //box constraint handling
        if( edge < DOFmin.get(dof) )
            edge = DOFmin.get(dof);
        else if( edge > DOFmax.get(dof) )
            edge = DOFmax.get(dof);


        //recursive attempts to get a good value when there are non-box constraints
        if(checkx){
            GCTol = 0.05;//very liberal
            boolean stillOOR = isOutOfRange(edge,dof);
            GCTol = 1e-10;//what it was before

            int maxNumRecursions = 10;

            if(stillOOR){
                if(numRecursions<10)//let's allow up to 10 recursions
                    edge = getEdgeDOFVal(edge,dof,checkx,numRecursions+1);
                else {
                    System.out.println("woah out of range edge value for CCDMinimizer, not fixed by "+maxNumRecursions+" recursions!!");

                    //getEdgeDOFVal(oorVal,dof);//INFINITE LOOP BUT WE'LL DO THIS FOR DEBUGGING!!

                    isOutOfRange(edge,dof);//For going into during debugging
                }
            }
        }

            //All GCs should be set up so that getNearestInRangeDOFVal returns the value at the constraints
            //closest to x.get(dof) (no farther than the starting value)
            //So, once we have moved edge to this boundary value, another GC
            //moving it closer to x.get(dof) cannot move it back out of range
            //so if we started with good x and
            //successively satisfy each constraint, the final value will be in range
            //(within numerical error)


        if(banZero)
            edge = enforceZeroBan(edge,dof);



        if(jumpOOR && checkx){//jumpOOR is intended for cases when we expect an in-range answer (i.e. checkx)
            double delta = jumpOORDelta;
            if(oorVal<edge)
                delta *= -1;

            for(int step=0; step<jumpOORNumSteps; step++){
                double DOFVal = oorVal + step*delta;
                if(!isOutOfRange(DOFVal,dof)){
                    //DOFVal is a legitimate possibility
                    if(objFcn.getValForDOF(dof,DOFVal)<objFcn.getValForDOF(dof, edge)){
                        edge = DOFVal;
                        break;
                    }
                }

                if(DOFVal<DOFmin.get(dof)||DOFVal>DOFmax.get(dof))
                    break;
                //we enforce box constraints throughout this process as usual...this is designed to jump
                //GC barriers not box constraint barriers
            }
        }


        return edge;
    }


    
    double enforceZeroBan(double edge, int dof){
        //Make sure x doesn't become 0 when edge is inserted at position dof
        if(edge==0){
            if(x.zDotProduct(x)==x.get(dof)*x.get(dof)){
                //If the DOF vector becomes 0 when this edge is put in

                if(x.get(dof)>0)
                    edge = 1e-10;//move the value up a little to avoid the 0, on the side nearer the original value
                else
                    edge = -1e-10;
            }
        }
        return edge;
    }


    
    public boolean compInitVals() {
       
        if(singleInitVal!=null){//use singleInitVal, w/o ellipse bounds
            x = getAnglesInRange(singleInitVal);

            for(int dof=0; dof<numDOFs; dof++){//return false if it's out of range
                if(isOutOfRange(x.get(dof),dof))
                    return false;
            }
            return true;
        }

        if(useCorners)
            compInitValsCorners();
        else
            compInitValsMiddle();

        return true;//indicates success
    }

    public void compInitValsMiddle() {
        //make initial values in the middle of the constraints
        double initVals[] = new double[numDOFs];
        for(int a=0; a<numDOFs; a++)
            initVals[a] = ( DOFmin.get(a) + DOFmax.get(a) ) / 2;

        x = DoubleFactory1D.dense.make(initVals);

        boolean badInit = satisfyNonBoxConstr();

        if(badInit){
            System.err.println("ERROR: Could not find feasible initial point for minimization.  GenCoord types:");
            for( GenCoord gc : nonBoxConstrGC )
                System.err.println(gc.type);
            new Exception().printStackTrace();
            System.exit(1);
        }
    }

    

    public void compInitValsCorners() {
        //make initial values at the corner with the lowest energy
        double firstCorner[] = new double[numDOFs];
        for(int a=0; a<numDOFs; a++)
            firstCorner[a] = DOFmin.get(a);
        
        x = DoubleFactory1D.dense.make(firstCorner);

        if(numDOFs == 0)//No corners to consider; not actually minimizing
            return;

        boolean upDown[] = new boolean[numDOFs];
        //indicates which DOFs are at their minimum (false) or maximum (true) value
        DoubleMatrix1D bestInitVals = null;

        double bestCornerE = Double.POSITIVE_INFINITY;
        boolean done = false;

        //Loop over corners
        while(!done){

            //Move in from each corner to the edge of the region allowed by the non-box constraints
            boolean badInit = satisfyNonBoxConstr();//badInit 

            //Use this corner's energy if it beats the others
            if( !badInit ){
                double curE = objFcn.getValue(x);
                if( curE < bestCornerE ){
                    bestCornerE = objFcn.getValue(x);
                    bestInitVals = x.copy();
                }
            }

            //Now go to the next corner
            int DOFToUpdate = numDOFs-1;

            while( upDown[DOFToUpdate] ){
                upDown[DOFToUpdate] = false;
                x.set(DOFToUpdate, DOFmin.get(DOFToUpdate));
                DOFToUpdate--;
                if(DOFToUpdate==-1){
                    done = true;
                    break;
                }
            }

            if(!done){
                upDown[DOFToUpdate] = true;
                x.set(DOFToUpdate,DOFmax.get(DOFToUpdate));
            }
        }

        if( Double.isInfinite(bestCornerE) ){
            System.err.println("ERROR: Could not find feasible initial point for minimization.  GenCoord types:");
            for( GenCoord gc : nonBoxConstrGC )
                System.err.println(gc.type);
            new Exception().printStackTrace();
            System.exit(1);
        }
        else{//Set the initial values to the best corner
            x = bestInitVals;
        }

    }

   

    DoubleMatrix1D getAnglesInRange(DoubleMatrix1D u){
        //Given a vector of this minimizer's DOFs, identify those that are angles
        //add multiples of 360 if possible to get them in the [DOFmin,DOFmax] range
        //angles are assumed to be bounded
        //either on both max and min side or not at all
        //this is for preparing initial values
        
        //we will move each out-of-range angle to the first equivalent angle after DOFmin
        DoubleMatrix1D ans = u.copy();
        
        for(int dof=0; dof<numDOFs; dof++){
            
            if( objFcn.isDOFAngle(dof) ){
            
                if(ans.get(dof)<DOFmin.get(dof)-0.001){
                    ans.set(dof,
                            DOFmin.get(dof) + (ans.get(dof)-DOFmin.get(dof))%(360) + 360 );
                }
                else if(ans.get(dof)>DOFmax.get(dof)+0.001){
                    ans.set(dof,
                            DOFmin.get(dof) + (ans.get(dof)-DOFmin.get(dof))%(360) );
                }
            }
        }

        return ans;
    }


    public void pickInitVals(ArrayList<DoubleMatrix1D> vecs){
        //from a list of possible initial-value vectors
        //that satisfy the constraints
        //pick the one lowest in energy and use it

        DoubleMatrix1D bestInitVals = null;
        double bestE = Double.POSITIVE_INFINITY;

        for(DoubleMatrix1D vec : vecs){

            double curE = objFcn.getValue(vec);
            if( curE < bestE ){
                bestE = curE;
                bestInitVals = vec;
            }
        }

        if( Double.isInfinite(bestE) ){
            System.err.println("ERROR: Could not find feasible initial point for minimization.  GenCoord types:");
            for( GenCoord gc : nonBoxConstrGC )
                System.err.println(gc.type);
            new Exception().printStackTrace();
            System.exit(1);
        }
        else{//Set the initial values to the best corner
            x = bestInitVals;
        }
    }


    private boolean satisfyNonBoxConstr() {
        //Move x in, dimension by dimension,
        //until we're in the region where the non-box constraints are satisfied
        //return true if they don't all get satisfied, e

        boolean badInit = false;

        //This procedure should at least work for common cases: no box constraints,
        //one linear-combination box constraint where the feasible region exists, etc.
        //Could maybe try repeating it one or more times for tougher cases?
        if(nonBoxConstrGC.length != 0){

            for(int dof=0;dof<numDOFs; dof++){
                double val = x.get(dof);
                if(isOutOfRange(val,dof)){
                    double edge = getEdgeDOFVal(val, dof);
                    //Since our point may not have been initially feasible, this step may not make it feasible
                    if( ! Double.isNaN(edge) )
                        x.set(dof, edge);
                }
            }

            //Now make sure the initial point is feasible
            for(int dof=0;dof<numDOFs; dof++){
                double val = x.get(dof);
                if(isOutOfRange(val,dof)){
                    badInit = true;
                    break;
                }
            }

        }

        return badInit;
    }

    

    public void compNonBoxConstrAffectingDOF(){
        //Compute nonBoxConstrAffectingDOF from nonBoxConstrIndices
        nonBoxConstrAffectingDOF = new int[numDOFs][];
        
        for(int dof=0; dof<numDOFs; dof++){

            ArrayList<Integer> affecting = new ArrayList<Integer>();
            for( int gcIndex=0; gcIndex<nonBoxConstrIndices.length; gcIndex++){
                for(int index : nonBoxConstrIndices[gcIndex]){
                    if(index == dof){
                        affecting.add(gcIndex);
                        break;
                    }
                }
            }

            nonBoxConstrAffectingDOF[dof] = new int[affecting.size()];
            for(int a=0; a<affecting.size(); a++)
                nonBoxConstrAffectingDOF[dof][a] = affecting.get(a);
        }
    }



    void rescaleValues() {
        double beta = rescalingGC.eval(x, rescalingIndices);

        int iter=0;

        while(beta==0){
            int randDOF = new java.util.Random().nextInt(numDOFs);
            double newVal = x.get(randDOF)+1e-5*Math.random();
            if(isOutOfRange(newVal,randDOF))
                newVal = getEdgeDOFVal(newVal,randDOF);
            x.set(randDOF,newVal);
            
            beta = rescalingGC.eval(x, rescalingIndices);
            iter++;

            if(iter>1e4)//5)
                return;
        }

        if(iter>0)
            System.out.println("Warning: rescalingGC=0 for CCD minimizer rescaling.  Trying to salvage with random perturbation.  "+iter+" iterations until beta!=0");

        if(beta==0)
            return;

        x.assign(Functions.mult(beta));

        //PROBLEM: this will result in -inf energy returned for beta==0...
    }



    double doRandMinCheck(double E){
        //This method is to check a minimum to make sure it's good to stop
        //We look at random vectors within an initial step size in each dimension
        //of x, and if any are lower in energy then we set x to that one and return its energy
        //otherwise we return the argument (the current energy of x)
        //and stop
        int numCheckPts = 100;
        int maxConsideredPts = 1000;//maximum number of points to consider
        double rangeFac = 10;//3;//1;//how many times step size to search over in each dimensions

        DoubleMatrix1D xCur = x.copy();

        int a=0;

        for(int c=0; c<maxConsideredPts && a<numCheckPts; c++){

            for(int dof=0; dof<numDOFs; dof++){
                double r = 2*Math.random()-1;//uniformly distributed from -1 to 1
                double change = rangeFac*objFcn.getInitStepSize(dof)*r;//getStepSize(dof,0)*r;
                x.set(dof, xCur.get(dof)+change);
            }


            //only consider in-range points
            boolean oor = false;//indicates point is out of range
            for(int dof=0; dof<numDOFs; dof++){
                if(isOutOfRange(x.get(dof),dof)){
                    oor = true;
                    break;
                }
            }

            if(!oor){
                double Ept = objFcn.getValue(x);

                if(Ept<E){//move x here
                    xCur = x.copy();
                    E = Ept;
                }

                a++;//getting numCheckPts actual samples
            }
        }

        x = xCur;
        objFcn.setDOFs(x);
        return E;
    }
    
    
    
    
    //Support for FBFB
    void recordCoordForFBFB(){
        //If using FBFB, we now have a favorable state from which we can linearly propagate
        //record coordinates for this
        if( objFcn instanceof ContSCObjFunction ){
            //this way lin approx will start at starting state
            //for(FullBBFlexBlock fbfb : ((ContSCObjFunction)objFcn).m.fbfb)
            //    fbfb.recordBlockActualCoords(false);
        }
    }
    
    
    void tryAltFBFBClosures(){
        //If there are FBFB's, try out any alternate solutions in their closures blocks
        if( objFcn instanceof ContSCObjFunction ){
            for( FullBBFlexBlock fbfb : ((ContSCObjFunction)objFcn).m.fbfb ){
                fbfb.tryAltClosures( ((ContSCObjFunction)objFcn).efunc );
            }
        }
    }

}