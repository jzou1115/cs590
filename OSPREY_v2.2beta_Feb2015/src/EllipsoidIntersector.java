
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.StrictlyConvexMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import java.util.ArrayList;
import java.io.Serializable;

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
//	EllipsoidIntersector.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class EllipsoidIntersector implements Serializable {
    //This is for finding intersections of ellipses
    //with each other and with box constraints
    //to use as initial values when we have EllipseBounds
    
    ArrayList<MinVolEllipse> ellipses;//these will be the outermost ellipses of all the EllipseBounds)
    //(since finite values of the EllipseBounds are the values inside these)
    //Because they come from the EllipseBounds, they are relative to the ellipse bound centers:
    ArrayList<DoubleMatrix1D> ellipseCenters;
    ArrayList<int[]> ellipseDOFs;//DOFNums for each ellipse
    
    DoubleMatrix1D boxConstraints[];
    CETEnergyFunction ef;
    
    double curInitPoint[] = null;
    
    int numFixedDOFs;//number of DOFs fixed by the box constraints to a single value
    //these will have to be taken out of minimization because they mess up the JOptimizer
    //the others are called free DOFs and here are mapping between them and the overall indices:
    int freeToOverallDOFs[] = null;
    int overallToFreeDOFs[] = null;
    
    
    int numDOFs;
    
    public EllipsoidIntersector(CETObjFunction of, DoubleMatrix1D ip){
        
        if(ip!=null)
            curInitPoint = ip.toArray();//initialize (perhaps to null)
        
        boxConstraints = of.constraints;
        ef = of.ef;
        ellipses = new ArrayList<MinVolEllipse>();
        ellipseCenters = new ArrayList<DoubleMatrix1D>();
        ellipseDOFs = new ArrayList<int[]>();
        
        for(ContETerm b : of.ef.terms){
            
            /*
            if(b instanceof EllipseBound){
                EllipseBound eb = (EllipseBound)b;
                
                int lastEllipseNum = eb.numEllipses-1;
                if(lastEllipseNum>-1){
                    while(eb.ellipses[lastEllipseNum]==null){
                        lastEllipseNum--;
                        if(lastEllipseNum==-1)//no ellipses: eb will always be 0, so don't need to account for it in feasible region
                            break;
                    }
                }
                if(lastEllipseNum>-1){
                    ellipses.add( eb.ellipses[lastEllipseNum] );
                    ellipseCenters.add(eb.center);
                    ellipseDOFs.add(eb.DOFNums);
                }
            }*/
            if(b.NCEllipse!=null){//our point must fall in the intersection of all non-null NCEllipses
                ellipses.add(b.NCEllipse);
                 ellipseCenters.add(b.center);//MAYBE KEEP AT 0
                 ellipseDOFs.add(b.DOFNums);
            }
        }
        
        numDOFs = boxConstraints[0].size();
        
        //count free DOFs
        numFixedDOFs = 0;
        overallToFreeDOFs = new int[numDOFs];
        int freeDOF = 0;
        for(int dof=0; dof<numDOFs; dof++){
            if(boxConstraints[1].get(dof)>boxConstraints[0].get(dof)){//free DOF
                overallToFreeDOFs[dof] = freeDOF;
                freeDOF++;
            }
            else{//fixed DOF
                overallToFreeDOFs[dof] = -1;
                numFixedDOFs++;
            }
        }
        
        freeToOverallDOFs = new int[numDOFs-numFixedDOFs];
        freeDOF = 0;
            
        for(int dof=0; dof<numDOFs; dof++){
            if(boxConstraints[1].get(dof)>boxConstraints[0].get(dof)){//free DOF
                freeToOverallDOFs[freeDOF] = dof;
                freeDOF++;
            }
        }
    }
    
    
    /*DoubleMatrix1D getFeasibleInitialPoint(){
        //Going to use JOptimizer to intersect the constraints
        //they are all convex so this is a convex optimization problem
        
        OptimizationRequest or = new OptimizationRequest();


        double[] q = new double[numDOFs];
        q[0] = 1;
        LinearMultivariateRealFunction objectiveFunction = 
                new LinearMultivariateRealFunction(q,0);//this is arbitrary
        or.setF0(objectiveFunction);
        
        
        //DEBUG: REMOVING ELLIPSOIDAL CONSTRAINTS
        
        
        int numConstr = 2*numDOFs;// + ellipses.size();
        ConvexMultivariateRealFunction[] inequalities 
                        = new ConvexMultivariateRealFunction[numConstr];
        
        int constrCount = 0;
        //linear constraints
        for(int dof=0; dof<numDOFs; dof++){
            double coeffs[] = new double[numDOFs];//coefficients for linear inequality for lower box constraint
            //(lower bound - x < 0)
            coeffs[dof] = -1;
            inequalities[constrCount] = new LinearMultivariateRealFunction(coeffs,boxConstraints[0].get(dof));
            constrCount++;
                        
            double coeffs2[] = new double[numDOFs];//upper box constraint
            coeffs2[dof] = 1;
            inequalities[constrCount] = new LinearMultivariateRealFunction(coeffs2,-boxConstraints[1].get(dof));
            constrCount++;
        }
        
        //ellipsoidal constraints
        for(int ell=0; ell<ellipses.size(); ell++){
            //DEBUG!!!
            //inequalities[constrCount] = new EllipsoidalConstraint(ellipses.get(ell),ellipseCenters.get(ell),ellipseDOFs.get(ell));
            constrCount++;
        }

        
        //DEBUG!!!
        double boxCenter[] = new double[numDOFs];
        for(int dof=0; dof<numDOFs; dof++)
            boxCenter[dof] = (boxConstraints[0].get(dof)+boxConstraints[1].get(dof))/2.;
        
        
        
        or.setFi(inequalities);

        //optimization problem
        //or.setNotFeasibleInitialPoint( new double[numParams] );
        or.setToleranceFeas(1.E-12);//as online
        or.setTolerance(1.E-12);
        //or.setInteriorPointMethod(JOptimizer.BARRIER_METHOD);

        //optimization
        JOptimizer opt = new JOptimizer();
        opt.setOptimizationRequest(or);
        try{
            int returnCode = opt.optimize();
            System.out.println("JOptimizer return code: "+returnCode);
        }
        catch(Exception e){
            System.err.print("JOptimizer error: "+e.getMessage());
            e.printStackTrace();//long stack traces in LS fits are causing trouble
            return null;//maybe this just means there's no initial point...
        }

        double v[] = opt.getOptimizationResponse().getSolution();
        return DoubleFactory1D.dense.make(v);
    }*/
    
    double[] removeFixedDOFs(double allDOFs[]){
        //convert from all DOFs to only free ones
        if(numFixedDOFs==0)
            return allDOFs;
        else {
            double[] ans = new double[numDOFs-numFixedDOFs];
            for(int dof=0; dof<numDOFs-numFixedDOFs; dof++)
                ans[dof] = allDOFs[freeToOverallDOFs[dof]];
            return ans;
        }
    }
    
    DoubleMatrix1D fillInFixedDOFs(DoubleMatrix1D freeDOFs){
        //reverse
        if(numFixedDOFs==0)
            return freeDOFs;
        else {
            DoubleMatrix1D ans = DoubleFactory1D.dense.make(numDOFs);
            for(int dof=0; dof<numDOFs; dof++) {
                int freeDOF = overallToFreeDOFs[dof];
                if(freeDOF==-1)
                    ans.set(dof,boxConstraints[0].get(dof));
                else
                    ans.set(dof,freeDOFs.get(freeDOF));
            }
            return ans;
        }
    }
    
    //version with only feasible initial points (add constraints one at a time)
    //since JOptimizer is apparently having trouble with getting initial points
    DoubleMatrix1D getFeasibleInitialPoint(){
        //Going to use JOptimizer to intersect the constraints
        //they are all convex so this is a convex optimization problem
        
        OptimizationRequest or = new OptimizationRequest();
        
        
        //start with just the box constraints, then minimize convexity constraints one by one
        //and, if we get a minimum below 0 (indicating feasibility), use that minimum
        //as the initial point for minimizing the next constraint...
        ConvexMultivariateRealFunction[] inequalities 
                        = new ConvexMultivariateRealFunction[2*(numDOFs-numFixedDOFs)];
        
        int constrCount = 0;
        //linear constraints
        int curFreeDOF = 0;//counter for non-fixed DOFs
        for(int dof=0; dof<numDOFs; dof++){
            
            if( boxConstraints[0].get(dof) < boxConstraints[1].get(dof) ){//only include free DOFs
                double coeffs[] = new double[numDOFs-numFixedDOFs];//coefficients for linear inequality for lower box constraint
                //(lower bound - x < 0)
                
                coeffs[curFreeDOF] = -1;
                inequalities[constrCount] = new LinearMultivariateRealFunction(coeffs,boxConstraints[0].get(dof));
                constrCount++;

                double coeffs2[] = new double[numDOFs-numFixedDOFs];//upper box constraint
                coeffs2[curFreeDOF] = 1;
                inequalities[constrCount] = new LinearMultivariateRealFunction(coeffs2,-boxConstraints[1].get(dof));
                constrCount++;
                
                curFreeDOF++;
            }
        }
        
        if(curInitPoint==null){//no possibly infeasible initial point given...start at box center
            curInitPoint = new double[numDOFs];
            //at first, the point is only required to be in the voxel
            //then we'll move it as needed to be in all the ellipses (one ellipse at a time)
            for(int dof=0; dof<numDOFs; dof++)
                curInitPoint[dof] = (boxConstraints[0].get(dof)+boxConstraints[1].get(dof))/2.;
        }
        else {
            //make sure our initial point is in the box
            //we'll then check it against ellipsoids & move it as needed
            for(int dof=0; dof<numDOFs; dof++){
                if( curInitPoint[dof] >= boxConstraints[1].get(dof) || curInitPoint[dof] <=  boxConstraints[0].get(dof) )
                    curInitPoint[dof] = (boxConstraints[0].get(dof)+boxConstraints[1].get(dof))/2.;
            }
        }
        
        curInitPoint = removeFixedDOFs(curInitPoint);
        
        or.setFi(inequalities);
        or.setInitialPoint(curInitPoint);
        //optimization problem
        
        /*or.setToleranceFeas(1.E-12);//as online
        or.setTolerance(1.E-12);*/
        //trying a little looser tolerance...might help with numerics when there are many ellipses, and 
        //anyway we just need to get into the feasible range
        or.setToleranceFeas(0.001);//chilling out...
        or.setTolerance(0.001);
        
        or.setInteriorPointMethod(JOptimizer.BARRIER_METHOD);

        
        
        for(int ell=0; ell<ellipses.size(); ell++){
            
            //constraint based on next ellipse (ellConstr<0 for inside of ellipse)
            StrictlyConvexMultivariateRealFunction ellConstr = new EllipsoidalConstraint(ellipses.get(ell),ellipseCenters.get(ell),ellipseDOFs.get(ell));
            
            if( ellConstr.value(curInitPoint) >= 0 ){//initial point currently on or outside this ellipse: need to move it inside
                
                or.setF0(ellConstr);

                //optimization
                JOptimizer opt = new JOptimizer();
                opt.setOptimizationRequest(or);
                try{
                    int returnCode = opt.optimize();
                    //System.out.println("JOptimizer return code: "+returnCode);
                    //kinda too much output...
                }
                catch(Exception e){
                    System.err.print("JOptimizer error: "+e.getMessage());
                    e.printStackTrace();
                    System.err.println("Outputting EllipsoidIntersector to ELdump.dat for inspection");
                    KSParser.outputObject(this, "ELdump.dat");
                    System.exit(1);//we can in principle handle this by just taking last level's FS term
                }

                curInitPoint = opt.getOptimizationResponse().getSolution();
                double bestConstrVal = ellConstr.value(curInitPoint);

                if(bestConstrVal>=0)//no feasible region
                    return null;

                or.setInitialPoint(curInitPoint);
            }
            
            if(ell==ellipses.size()-1)//last ellipse: now curInitPoint is inside all the ellipses
                return fillInFixedDOFs(DoubleFactory1D.dense.make(curInitPoint));

            
            //add the constraint for this ellipsoid to the optimization problem
            ConvexMultivariateRealFunction[] newInequalities 
                        = new ConvexMultivariateRealFunction[inequalities.length+1];
            System.arraycopy(inequalities,0,newInequalities,0,inequalities.length);
            newInequalities[inequalities.length] = ellConstr;
            inequalities = newInequalities;
            or.setFi(inequalities);
        }
        
        //shouldn't get here!!!
        return null;
        
    }
    
    
    
    
    
    private class EllipsoidalConstraint implements StrictlyConvexMultivariateRealFunction {
        //this is a constraint demanding that we be within the ellipsoid
        MinVolEllipse ellipse;
        double[][] hess;//The Hessian is always the ellipse's matrix, mapped into the right DOF components
        DoubleMatrix1D center;
        int DOFNums[];

        public EllipsoidalConstraint(MinVolEllipse ellipse, DoubleMatrix1D center, int[] DOFNums) {
            this.ellipse = ellipse;
            this.center = center;
            this.DOFNums = DOFNums;
            
            //Hessian for the involved DOFs only
            DoubleMatrix2D hessRed = ellipse.A.copy();
            hessRed.assign(Functions.mult(2));
            
            //Map hess into the right components of the overall DOF vector space
            //the other components of the full-DOF Hessian are 0
            hess = new double[numDOFs-numFixedDOFs][numDOFs-numFixedDOFs];
            for(int dof=0; dof<DOFNums.length; dof++){
                for(int dof2=0; dof2<DOFNums.length; dof2++){
                    int freeDOF1 = ef.revDOFNums.get(DOFNums[dof]);
                    int freeDOF2 = ef.revDOFNums.get(DOFNums[dof2]);
                    
                    if(numFixedDOFs>0){//there are fixedDOFs
                        freeDOF1 = overallToFreeDOFs[freeDOF1];
                        if(freeDOF1!=-1){//need this row
                            freeDOF2 = overallToFreeDOFs[freeDOF2];
                            if(freeDOF2!=-1)
                                hess[freeDOF1][freeDOF2] = hessRed.get(dof,dof2);
                        }
                    }
                    else
                        hess[freeDOF1][freeDOF2] = hessRed.get(dof,dof2);
                }
            }
        }
        
        
        @Override
        public int getDim(){
            return numDOFs-numFixedDOFs;
        }
        
        
        @Override
        public double value(double[] x){
            DoubleMatrix1D z = getEllipseInput(x);
            return ellipse.getScaledDistFromCenter(z) - 1;
        }
        
        
        @Override
        public double[] gradient(double[] x){
            DoubleMatrix1D z = getEllipseInput(x);
            z.assign(ellipse.nc,Functions.plus);
            z = Algebra.DEFAULT.mult(ellipse.A, z);
            z.assign(Functions.mult(2));
            
            double ans[] = new double[numDOFs - numFixedDOFs];
            for(int dof=0; dof<DOFNums.length; dof++){
                
                int freeDOF = ef.revDOFNums.get(DOFNums[dof]);
                if(numFixedDOFs>0){//there are fixedDOFs
                    freeDOF = overallToFreeDOFs[freeDOF];
                    if(freeDOF!=-1)
                        ans[freeDOF] = z.get(dof);
                }
                else
                    ans[freeDOF] = z.get(dof);
            }
            
            
            
            //TEST!!!!!!
            /*double gradStep = 1e-4;
            double numGrad[] = new double[numDOFs];
            for(int dof=0; dof<numDOFs; dof++){
                double newVec[] = x.clone();
                newVec[dof] += gradStep;
                
                double newVecn[] = x.clone();
                newVecn[dof] -= gradStep;
                numGrad[dof] = (value(newVec)-value(newVecn))/(2*gradStep);
            }*/
            
            
            
            return ans;
        }
        
        
        DoubleMatrix1D getEllipseInput(double[] x){
            //convert x to the input form for this particular ellipse
            DoubleMatrix1D overallDOFVals = fillInFixedDOFs(DoubleFactory1D.dense.make(x));
            DoubleMatrix1D z = ef.getTermDOFVals(overallDOFVals, DOFNums, null);
            z.assign(center,Functions.minus);//z-center
            return z;
        }
        
        
        @Override
        public double[][] hessian(double[] x){    
            
            
                       //TEST!!!!!!
            /*double gradStep = 1e-4;
            double numHess[][] = new double[numDOFs][numDOFs];
            for(int dof=0; dof<numDOFs; dof++){
                for(int dof2=0; dof2<numDOFs; dof2++){
                    
                    double newVecpp[] = x.clone();
                    newVecpp[dof] += gradStep;
                    newVecpp[dof2] += gradStep;
                    
                    double newVecpn[] = x.clone();
                    newVecpn[dof] += gradStep;
                    newVecpn[dof2] -= gradStep;
                    
                    double newVecnp[] = x.clone();
                    newVecnp[dof] -= gradStep;
                    newVecnp[dof2] += gradStep;
                    
                    double newVecnn[] = x.clone();
                    newVecnn[dof] -= gradStep;
                    newVecnn[dof2] -= gradStep;
                    
                    numHess[dof][dof2] = (value(newVecpp)-value(newVecpn)-value(newVecnp)+value(newVecnn))/(4*gradStep*gradStep);
                }
            } 
            
            
            double numHess2[][] = new double[numDOFs][numDOFs];
            for(int dof=0; dof<numDOFs; dof++){
                
                double newVec[] = x.clone();
                newVec[dof] += gradStep;
                double gradp[] = gradient(newVec);
                
                double newVecn[] = x.clone();
                newVecn[dof] -= gradStep;
                double gradn[] = gradient(newVecn);
                
                for(int dof2=0; dof2<numDOFs; dof2++)
                    numHess2[dof][dof2] = (gradp[dof2]-gradn[dof2])/(2*gradStep);
            } */
            
            
            return hess;
        }
    }
    
    
    
}
