
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
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
//	FullBBFlexBlock.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class FullBBFlexBlock implements Serializable {
    //this vicious monster is capable of moving the backbone in any direction!
    //It can grab the backbone at arbitrary points in the middle and pull in the desired direction
    //thus allowing comprehensive search of backbone space!!!!
    
    //we assume the backbone's Atom.coord are the starting coordinates for these backbone adjustments
    
    Molecule m;//the unsuspecting molecule, which will soon suffer changes in its very backbone!
    int startRes;
    int endRes;
    //the section from CA of startRes to CA of endRes is flexible
    
    
    int freePepPlanes[];//These numbers are molecule residue numbers for the first CA of each peptide plane
        
    int closureBlock[];//Molecule residue numbers i such that Cai is the middle of a block needing closure
    //using the new octic algorithm
    
    LoopClosureAdjustment blockLCAs[];//LCAs for closureBlocks.  Stored only so we can go through different
    //solutions if needed
    
    int blockNumSolns[];//number of solutions for each closure block
    
            
    double[] genChi1;//(Generalized) chi1 angles for all the residues in this.  We'll record them
    //before we apply new DOF values and then set them the same once we get a BB conf in place
    //(so that sidechain dihedral motions are kept separate from BB motions, as in Perturbation)
    
    
    BBFreeDOF freeDOFs[];
    //possible types:
    static int ANCHORDIH=0, FIRSTCAADJ=1, SECONDCAADJ=2, PEPPLANESPIN=3;
    
    int numFreeDOFs;
    
    double[][][] blockActualCoordOld, blockActualCoordNew;//for recording mainchain atom
    //actualCoordinates in closure blocks

    double maxCATrans = Double.POSITIVE_INFINITY;//when we do loop closures,
    //we'll make sure the new CA coordinates are no further than this from starting points (in each dim)
    //it makes sense for this to be the same as the bounds on the FIRSTCAADJ translation parameters

    public FullBBFlexBlock(Molecule m, int startRes, int endRes, boolean anchorPreference, double maxCATrans) {
        //if we need anchor dihedrals on one anchor, anchorPreference=true chooses starting anchor
        
        this.m = m;
        this.startRes = startRes;
        this.endRes = endRes;
        this.maxCATrans = maxCATrans;
        
        int totNumPepPlanes = endRes-startRes;
        int numFreePlanes = (endRes-startRes-2)/3;
        
        freePepPlanes = new int[numFreePlanes];
        
        closureBlock = new int[numFreePlanes+1];//must have 1 of these before and after each free plane
                
        if(numFreePlanes<1)
            throw new RuntimeException("ERROR: FullBBFlexBlock too small.  startRes: "+startRes
                    +" endRes: "+endRes);
        
        
        freePepPlanes[0] = startRes+2;
        if( totNumPepPlanes%3==1 || (anchorPreference &&  totNumPepPlanes%3==0) )
            freePepPlanes[0]++;//anchor dihedrals on first anchor
        
        closureBlock[0] = freePepPlanes[0]-1;
        
        //OK now we have free pep planes and closure blocks repeating every 3 residues
        for(int pl=1; pl<numFreePlanes; pl++)
            freePepPlanes[pl] = freePepPlanes[pl-1]+3;
        
        for(int bl=1; bl<=numFreePlanes; bl++)
            closureBlock[bl] = closureBlock[bl-1]+3;
        
        
        numFreeDOFs = 6*numFreePlanes;//actual free pep plane DOFs
        if(totNumPepPlanes%3==1)//two anchor dihedrals
            numFreeDOFs += 4;
        else if(totNumPepPlanes%3==0)//one anchor dihedral
            numFreeDOFs += 2;
        
        freeDOFs = new BBFreeDOF[numFreeDOFs];
        for(int pl=0; pl<numFreePlanes; pl++){//free pep plane DOFs
            
            //translations
            for(int dim=0; dim<3; dim++)
                freeDOFs[6*pl+dim] = new BBFreeDOF(FIRSTCAADJ,0,dim,pl);
            
            //rotations
            freeDOFs[6*pl+3] = new BBFreeDOF(SECONDCAADJ,0,0,pl);
            freeDOFs[6*pl+4] = new BBFreeDOF(SECONDCAADJ,0,1,pl);
            freeDOFs[6*pl+5] = new BBFreeDOF(PEPPLANESPIN,0,0,pl);
        }
        
        //make freeDOFs for anchor dihedrals
        int count = 6*numFreePlanes;
        
        if( totNumPepPlanes%3==1 || (anchorPreference &&  totNumPepPlanes%3==0) ){
            //anchor dihedrals on first anchor
            freeDOFs[count] = new BBFreeDOF(ANCHORDIH,0,0,0);
            count++;
            freeDOFs[count] = new BBFreeDOF(ANCHORDIH,1,0,0);
            count++;
        }
        
        if( totNumPepPlanes%3==1 || ( (!anchorPreference) &&  totNumPepPlanes%3==0) ){
            //anchor dihedrals on 2nd anchor
            freeDOFs[count] = new BBFreeDOF(ANCHORDIH,2,0,0);
            count++;
            freeDOFs[count] = new BBFreeDOF(ANCHORDIH,3,0,0);
            count++;
        }
        
        
        
        blockLCAs = new LoopClosureAdjustment[numFreePlanes+1];
        blockNumSolns = new int[numFreePlanes+1];
        blockActualCoordOld = new double[numFreePlanes+1][9][3];//record 9 atoms per block
        blockActualCoordNew = new double[numFreePlanes+1][9][3];
        
        recordBlockActualCoords(false);//record coordinates so we can do linear propagation
        
        genChi1 = new double[endRes-startRes+1];
    }
    
    
  
    

    
    
    
    protected class BBFreeDOF implements Serializable {
        int type;//one of the 4 types
        int anchorDihNum;//for type ANCHORDIH: 0 and 1 for beginning of block, 2 and 3 for end (in sequence order)
        int coordNum;//For CA adjustment: 
        int freePlaneNum;//which free peptide plane this free DOF is for

        public BBFreeDOF(int type, int anchorDihNum, int coordNum, int freePlaneNum) {
            this.type = type;
            this.anchorDihNum = anchorDihNum;
            this.coordNum = coordNum;
            this.freePlaneNum = freePlaneNum;
        }
    }
    
    
    //we assume that anchordih's are applied starting at the anchors (so #1 after #0 and #2 after #3)
    //we also assume that within each free peptide plane, we apply the first CA adjustments 
    //(translating the whole plane starting from the Atom.coord)
    //and then the second CA adjustments
    //and then the pep plane spin
    //these orderings are set when we create this FullBBFlexBlock
    
    
    boolean setFreeDOFs(DoubleMatrix1D x){
        //set all the free DOFs
        //return true normally but false if no way to close loop
        
        for(int res=startRes; res<=endRes; res++)
            genChi1[res-startRes] = Perturbation.getGenChi1(m, res);
        
        backupCoords();
        
        //recordBlockActualCoords(false);
        //record actualCoordinates 
        //will use for linear propagation
        //we will want to do linear propagating only starting from steps known to be favorable
        //so this will be done outside
        
        
        for(int dof=0; dof<x.size(); dof++){
            
            BBFreeDOF curDOF = freeDOFs[dof];
            double val = x.get(dof);
            
            if(curDOF.type==ANCHORDIH)
                applyAnchorDih(curDOF.anchorDihNum,val);
            else if (curDOF.type==FIRSTCAADJ)
                translatePepPlane(curDOF.freePlaneNum,curDOF.coordNum,val);
            else if (curDOF.type==SECONDCAADJ)
                rotatePepPlane(curDOF.freePlaneNum,curDOF.coordNum,val);
            else//PEPPLANESPIN
                rotatePepPlane(curDOF.freePlaneNum,2,val);
        }
        
        
        recordBlockActualCoords(true);
            
        boolean success = true;//able to close chain
        
        for(int block=0; block<closureBlock.length; block++){
            success = success && calcSolns(block);//get the possible solutions...NOTE PROBABLY FASTEST IF DONT USE ALL OF THESE
            
            if(!success){//no way to close this block...put back the actualCoordinates and return false
                revertCoords();
                return false;
            }
            
            int bestSoln = pickSoln(block);//get the closest to the linearly propagated value
            
            if(bestSoln==-1){//no way to close the block with CAMaxTrans filter satisfied
                revertCoords();
                return false;
            }
            
            applySoln(block,bestSoln);//move the backbones (linked to sidechains as rigid body)
        }
        
        //now the backbones are in the correct conformations
        //get all the sidechains into the correct conformations too
        for(int res=startRes; res<=endRes; res++){
            if(!fixSC(res)){//conformation isn't OK, because we can't close some sidechain(s)
                //note we still declare failure even if there are other backbone options,
                //because this is the backbone option that follows continuously in our minimization
                revertCoords();
                return false;
            }
        }
        
        return true;
    }
    
    
    
    
    //backup mechanism in case we can't close the loop
    double[][] backupCoord;
    
    void backupCoords(){
        backupCoord = new double[endRes-startRes+1][];
        for(int res=0; res<backupCoord.length; res++)
            backupCoord[res] = storeResActualCoord(startRes+res);
    }
    
    double[] storeResActualCoord(int molResNum){
        Residue curRes = m.residue[molResNum];
        double ans[] = new double[3*curRes.numberOfAtoms];
        System.arraycopy( m.actualCoordinates, 3*curRes.atom[0].moleculeAtomNumber, ans, 0,
                3*curRes.numberOfAtoms );
        return ans;
    }
    
    
    void revertCoords(){
        for(int res=0; res<backupCoord.length; res++)
            restoreResActualCoord(startRes+res,backupCoord[res]);
    }
    
    void restoreResActualCoord(int molResNum, double[] coords){
        Residue curRes = m.residue[molResNum];
        System.arraycopy(coords, 0, m.actualCoordinates, 3*curRes.atom[0].moleculeAtomNumber,
                3*curRes.numberOfAtoms );
    }
    
    
    
    //application of free DOFs
    //we move around sidechains as rigid bodies; then, as long at the CA-CB distance and CA-CB-X angles
    //are kept as in the original structure, Molecule.idealizeSidechain can be used to get the correct geometry
    //(and we will need to set genChi1 too of course)
    //Also, we move atoms in the closure blocks to prepare them for the rotations in LoopClosureAdjustment.applyTripeptideRot
    //that is, the middle carbonyl and everything beyond it stays connected to the peptide plane after the closure block
    //and the rest stays connected to the plane before the closure block
    //this way, applyTripeptideRot will reconnect the chain and restore correct geometry
    
    void applyAnchorDih(int anchorDihNum, double val){
        //we aim to readjust either the first or the last peptide plane of this FullBBFlexBlock
        
        if(anchorDihNum==0){//first dihedral building on starting anchor
            int C0 = m.residue[startRes-1].getAtomNameToMolnum("C");
            int N1 = m.residue[startRes].getAtomNameToMolnum("N");
            int CA1 = m.residue[startRes].getAtomNameToMolnum("CA");
            int C1 = m.residue[startRes].getAtomNameToMolnum("C");
            
            //we rotate the section up to the middle carbonyl of the ensuing closure block
            int atomList1[] = m.residue[startRes+1].getAtomList(true, true, true, true);
            int atomList2[] = m.residue[startRes+2].getAtomList(true, true, true, false);
            int ONum = m.residue[startRes].getAtomNameToMolnum("O");//we also need to move carbonyl oxygen of startRes
            
            int atomList[] = LoopClosureAdjustment.concatenateArrays( new int[] {ONum}, atomList1, atomList2 );
            m.setTorsion(C0, N1, CA1, C1, val, atomList, atomList.length, false);
        }
        else if(anchorDihNum==1){//second
            int N1 = m.residue[startRes].getAtomNameToMolnum("N");
            int CA1 = m.residue[startRes].getAtomNameToMolnum("CA");
            int C1 = m.residue[startRes].getAtomNameToMolnum("C");
            int N2 = m.residue[startRes+1].getAtomNameToMolnum("N");
            
            int atomList1[] = m.residue[startRes+1].getAtomList(false, true, true, true);
            int atomList2[] = m.residue[startRes+2].getAtomList(true, true, true, false);
            int HNum = m.residue[startRes+1].getAtomNameToMolnum("H");//amide is accounted for separately
            
            int atomList[] = LoopClosureAdjustment.concatenateArrays( new int[] {HNum}, atomList1, atomList2 );
            m.setTorsion(N1, CA1, C1, N2, val, atomList, atomList.length, false);
        }
        else if(anchorDihNum==3){//first dihedral on ending anchor
            int N2 = m.residue[endRes+1].getAtomNameToMolnum("N");
            int C1 = m.residue[endRes].getAtomNameToMolnum("C");
            int CA1 = m.residue[endRes].getAtomNameToMolnum("CA");
            int N1 = m.residue[endRes].getAtomNameToMolnum("N");
            
            int atomList1[] = m.residue[endRes-1].getAtomList(true, true, true, true);
            int atomList0[] = m.residue[endRes-2].getAtomList(false, false, false, true);//we just move the carbonyl of the closure block's middle residue
            int HNum = m.residue[endRes].getAtomNameToMolnum("H");//amide is accounted for separately
            
            int atomList[] = LoopClosureAdjustment.concatenateArrays( atomList0, atomList1, new int[] {HNum} );

            m.setTorsion(N2, C1, CA1, N1, val, atomList, atomList.length, false);
        }
        else if(anchorDihNum==2){//second
            int C1 = m.residue[endRes].getAtomNameToMolnum("C");
            int CA1 = m.residue[endRes].getAtomNameToMolnum("CA");
            int N1 = m.residue[endRes].getAtomNameToMolnum("N");
            int C0 = m.residue[endRes-1].getAtomNameToMolnum("C");
            
            int atomList1[] = m.residue[endRes-1].getAtomList(true, true, true, false);
            int atomList0[] = m.residue[endRes-2].getAtomList(false, false, false, true);
            int ONum = m.residue[endRes-1].getAtomNameToMolnum("O");//carbonyl is accounted for separately
            
            int atomList[] = LoopClosureAdjustment.concatenateArrays( atomList0, atomList1, new int[] {ONum} );
            m.setTorsion(C1, CA1, N1, C0, val, atomList, atomList.length, false);
        }
    }
    
    
    double measureAnchorDih(int anchorDihNum, boolean useActualCoord){
        //same deal but measuring current dihedrals (or Atom.coord dihedrals, if useActualCoord is false)
        
        if(anchorDihNum==0){//first dihedral building on starting anchor
            int C0 = m.residue[startRes-1].getAtomNameToMolnum("C");
            int N1 = m.residue[startRes].getAtomNameToMolnum("N");
            int CA1 = m.residue[startRes].getAtomNameToMolnum("CA");
            int C1 = m.residue[startRes].getAtomNameToMolnum("C");
            
            return measureDihedral(C0, N1, CA1, C1, useActualCoord);
        }
        else if(anchorDihNum==1){//second
            int N1 = m.residue[startRes].getAtomNameToMolnum("N");
            int CA1 = m.residue[startRes].getAtomNameToMolnum("CA");
            int C1 = m.residue[startRes].getAtomNameToMolnum("C");
            int N2 = m.residue[startRes+1].getAtomNameToMolnum("N");
            
            return measureDihedral(N1, CA1, C1, N2, useActualCoord);
        }
        else if(anchorDihNum==3){//first dihedral on ending anchor
            int N2 = m.residue[endRes+1].getAtomNameToMolnum("N");
            int C1 = m.residue[endRes].getAtomNameToMolnum("C");
            int CA1 = m.residue[endRes].getAtomNameToMolnum("CA");
            int N1 = m.residue[endRes].getAtomNameToMolnum("N");
            
            return measureDihedral(N2, C1, CA1, N1, useActualCoord);
        }
        else if(anchorDihNum==2){//second
            int C1 = m.residue[endRes].getAtomNameToMolnum("C");
            int CA1 = m.residue[endRes].getAtomNameToMolnum("CA");
            int N1 = m.residue[endRes].getAtomNameToMolnum("N");
            int C0 = m.residue[endRes-1].getAtomNameToMolnum("C");
            
            return measureDihedral(C1, CA1, N1, C0, useActualCoord);
        }
        else
            throw new RuntimeException("ERROR: Unrecognized anchor dihedral #: "+anchorDihNum);
    }
    
    
    double measureDihedral(int at1, int at2, int at3, int at4, boolean useActualCoord){
        //given 4 molecule atoms numbers, return their dihedral,
        //for actualCoordinates or Atom.coord as specified
        if(useActualCoord)
            return m.getTorsion(at1,at2,at3,at4);
        else
            return m.atom[at4].torsion( m.atom[at1], m.atom[at2], m.atom[at3] );
    }
    
    
    //note: we can likely speed up {translate,rotate}PepPlane
    //by applying a single 3-D translation and then a single rotation matrix
    //this may not be a bottleneck though
    
    
    void translatePepPlane(int freePlaneNum,int coordNum,double val){
        //translate the given free peptide plane along the given coordinate by the given amount
        int atomList[] = getFreePepTerritory(freePlaneNum);
        double trVec[] = new double[3];
        
        double relVal = val - measurePepPlaneTrans(freePlaneNum,coordNum);
        
        trVec[coordNum] = relVal;
        m.translateAtomList(atomList, trVec, false, false);
    }
    
    
    double measurePepPlaneTrans(int freePlaneNum,int coordNum){
        //measure current translation
        return m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]].getAtomNameToMolnum("CA"))[coordNum]
                - m.residue[freePepPlanes[freePlaneNum]].getAtomByName("CA").coord[coordNum];
    }
    
    
    
    void rotatePepPlane(int freePlaneNum, int coordNum, double val){
        //same deal but rotating
        //ranges of angles are 0 to 180, -180 to 180, 0 to 360
        
        int atomList[] = getFreePepTerritory(freePlaneNum);
        
        
        RotMatrix r = new RotMatrix();
        
        double axis[];
        
        //all rotations are about the first CA
        double actualCA1[] = m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]].getAtomNameToMolnum("CA"));
        double actualCA2[] = m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]+1].getAtomNameToMolnum("CA"));
        double[] actualCACA = r.subtract(actualCA2,actualCA1);

        Atom CA1 = m.residue[freePepPlanes[freePlaneNum]].getAtomByName("CA");
        Atom CA2 = m.residue[freePepPlanes[freePlaneNum]+1].getAtomByName("CA");
        double CACAVec[] = r.subtract( CA2.coord, CA1.coord );
        Atom N2 = m.residue[freePepPlanes[freePlaneNum]+1].getAtomByName("N");
        double CANVec[] = r.subtract( N2.coord, CA1.coord );

        
        //m.rotateAtomList will normalize the axis
        if(coordNum<2){//getting second CA in place.
            //rotation axes are based on the original peptide plane normal
                      

            if(coordNum==0){//rotate around normal.  If we let CA1 be origin and CA2 start on polar axis, this sets latitude...
                double curLatitude = 180*r.getAngle(CACAVec,actualCACA)/Math.PI;
                val -= curLatitude;//only need to rotate enough to get to the right latitude from curLatitude
                
                axis =  r.cross(CACAVec,actualCACA);//move directly towards or away from pole
                
                if(r.normsq(axis)==0)//if at pole, move toward original nitrogen position
                    axis = r.cross(CACAVec, CANVec);
            }
            else{//...and this sets longitude
                axis = CACAVec;
                double curLongitude = 180*r.getLongitude(actualCACA,CACAVec,CANVec)/Math.PI;
                val -= curLongitude;
            }
        }
        else {//for coordNum==2, we already have second CA in place
            //so we just spin around the CA1-CA2 axis (in the actualCoordinates)
            axis = actualCACA;
            
            double actualN2[] = m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]+1].getAtomNameToMolnum("N"));
            double[] actualCAN = r.subtract(actualN2,actualCA1);
            
            double curAngle;
            
            if( r.normsq(r.subtract(actualCACA, CACAVec))<1e-12 ){
                curAngle = 180*r.getLongitude(actualCAN,actualCACA,CANVec)/Math.PI;//then we can just use CANVec for longitude!
                if(curAngle<=0)//get into right range
                    curAngle += 360;
            }
            else
                curAngle = 180 + 180*r.getLongitude(actualCAN,actualCACA,CACAVec)/Math.PI;
            //we'll make the range of this (0,360]
            //we use the fact that for spin angle=0, original (only translated) CA2
            //is in same plane as actual N, CA2, and CA1, so CACAVec is at longitude 180
                
            
            val -= curAngle;
        }
        
        m.rotateAtomList( atomList, axis[0], axis[1], axis[2], actualCA1[0], actualCA1[1], 
                actualCA1[2], val, false );
    }
    
    
    
    double measurePepPlaneRot(int freePlaneNum, int coordNum){
        //structured just like rotatePepPlane, but measuring current values
        
        RotMatrix r = new RotMatrix();
                
        //all rotations are about the first CA
        double actualCA1[] = m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]].getAtomNameToMolnum("CA"));
        double actualCA2[] = m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]+1].getAtomNameToMolnum("CA"));
        double[] actualCACA = r.subtract(actualCA2,actualCA1);

        Atom CA1 = m.residue[freePepPlanes[freePlaneNum]].getAtomByName("CA");
        Atom CA2 = m.residue[freePepPlanes[freePlaneNum]+1].getAtomByName("CA");
        double CACAVec[] = r.subtract( CA2.coord, CA1.coord );
        Atom N2 = m.residue[freePepPlanes[freePlaneNum]+1].getAtomByName("N");
        double CANVec[] = r.subtract( N2.coord, CA1.coord );

        
        if(coordNum==0)
            return 180*r.getAngle(CACAVec,actualCACA)/Math.PI;
        else if (coordNum==1)
            return 180*r.getLongitude(actualCACA,CACAVec,CANVec)/Math.PI;
        else {//coordNum==2
            //so we just spin around the CA1-CA2 axis (in the actualCoordinates)

            
            double actualN2[] = m.getActualCoord(m.residue[freePepPlanes[freePlaneNum]+1].getAtomNameToMolnum("N"));
            double[] actualCAN = r.subtract(actualN2,actualCA1);
            
            double curAngle;
            
            if( r.normsq(r.subtract(actualCACA, CACAVec))<1e-12 ){
                curAngle = 180*r.getLongitude(actualCAN,actualCACA,CANVec)/Math.PI;//then we can just use CANVec for longitude!
                if(curAngle<=0)//get into right range
                    curAngle += 360;
            }
            else
                curAngle = 180 + 180*r.getLongitude(actualCAN,actualCACA,CACAVec)/Math.PI;

            return curAngle;
        }
    }
    
    
    
    
    int[] getFreePepTerritory(int freePlaneNum){
        //List of atoms that move along with the given free peptide plane
        //those that aren't part of the plane (they're part of closure blocks instead)
        //will then be rotated about one of the plane's CA's
        //by LoopClosureAdjustment.applyTripeptideRot
        //to reconnect the mainchains
        
        //start with the carbonyl before the peptide plane begins...  (middle res of preceding closure block)
        int atomList0[] = m.residue[freePepPlanes[freePlaneNum]-1].getAtomList(false, false, false, true);
        int atomList1[] = m.residue[freePepPlanes[freePlaneNum]].getAtomList(true, true, true, true);
        int atomList2[] = m.residue[freePepPlanes[freePlaneNum]+1].getAtomList(true, true, true, true);
        //middle residue of next closure block: move everything but carbonyl
        int atomList3[] = m.residue[freePepPlanes[freePlaneNum]+2].getAtomList(true, true, true, false);    
        
        return LoopClosureAdjustment.concatenateArrays(atomList0,atomList1,atomList2,atomList3);
    }
    
    
    


    //(lets say we can set psi/phi for second CA of freePepPlane, rel to first, start on x-axis)


    DoubleMatrix1D recordFreeDOFs(){
        //so like to recapitulate the L54W motion we can deepcopy m (from wt system),
        //copy the desired mut peptide planes from W to actualCoordinates, and measure
        //can be nice for loop grafting w/ minimal change
        
        //ALSO USE TO CHECK THAT DOFs ARE APPLIED PROPERLY
        
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(freeDOFs.length);
        
        for(int dof=0; dof<freeDOFs.length; dof++){
            
            BBFreeDOF curDOF = freeDOFs[dof];
            double val;
            
            if(curDOF.type==ANCHORDIH)
                val = measureAnchorDih(curDOF.anchorDihNum, true);
            else if (curDOF.type==FIRSTCAADJ)
                val = measurePepPlaneTrans(curDOF.freePlaneNum,curDOF.coordNum);
            else if (curDOF.type==SECONDCAADJ)
                val = measurePepPlaneRot(curDOF.freePlaneNum,curDOF.coordNum);
            else
                val = measurePepPlaneRot(curDOF.freePlaneNum,2);
            
            ans.set(dof,val);
        }
        
        return ans;
    }



    
    boolean calcSolns(int blockNum){
        //Based on the anchor atom coords for the specified closure block currently specified
        //in the actualCoordinates, get the solutions
        //and put them in blockSolns
        //THIS INITIAL VERSION JUST USES CSJD
        //WE CAN LATER INTRODUCE VOXEL-BASED BRACKETING OF POLYNOMIAL ROOT
        //AND EVEN EPIC FITTING OF SOLN (in some appropriate form)
        //WE'll USE A LOOPCLOSUREADJUSTMENT FOR CONVENIENCE EVEN THOUGH IT ISN'T CALLED LIKE A USUAL PERTURBATION
        //IT WILL BE STORED JUST WHILE WE HAVE THE CURRENT POSE OF THE PEPTIDE PLANES
        //ONE EITHER SIDE OF THE CLOSURE BLOCK
        
        int resList[] = new int[] {closureBlock[blockNum]-1, closureBlock[blockNum], closureBlock[blockNum]+1};
        blockLCAs[blockNum] = new LoopClosureAdjustment( "LOOP CLOSURE ADJUSTMENT", m, resList );
        //the difference from a usual LCA is we want the original (Atom.coord) bond lengths, angles, and omegas...
        //the chain may be broken in the actualCoordinates
        //NOTE: this would also be a good place to insert root bounds for the TripeptideClosure (based on middle-CA voxel, etc.)
        blockLCAs[blockNum].tc = new TripeptideClosure(m,resList,false);
        blockLCAs[blockNum].predecessors = new int[0];
        blockLCAs[blockNum].updateTC = false;//want to keep the Atom.coord-based tc
        //otherwise we'll end up with the weird bond lengths and angles that come from un-closed free pep plane motions!
        blockNumSolns[blockNum] = blockLCAs[blockNum].calcSolns();
        
        return (blockNumSolns[blockNum]>0);//return whether we found any solutions or not
        
        //if considering all solns, may also want to screen for good voxels?
    }
    
    
    int pickSoln(int block){
        //figure out which solution for the given closure block is most appropriate for the
        //associated change in free DOFs
        //meaning, which is closest to the linear approximation in a least-squares sense
        //for finite maxCATrans, we also filter out solutions that move the middle CA much
        //relative to initial (Atom.coord) position
        //If everything is filtered out, return -1
        
        double[][] linApprox = getBlockLinApprox(block);
        
        double minSSD=Double.POSITIVE_INFINITY;//maximum sum of square deviations
        int ans = -1;
        
        for(int soln=0; soln<blockNumSolns[block]; soln++){
            double ssd = getBlockLinApproxDeviation(block,soln,linApprox);//will be inf if CAMaxTrans filter violated!
            if(ssd<minSSD){
                minSSD = ssd;
                ans = soln;
            }
        }
            
        return ans;
    }
    
    
    double[][] getBlockLinApprox(int block){
        //linearly approximate the atom coords in the block, as given in blockAtoms
        //this is based on taking the linear relationship between time derivatives
        //and using this with the finite motion of the free atoms around the block
        
        //build the linear system
        //variables will be the coords for atoms in blockAtoms, one after the other
        //there are 9 atoms and thus 27 variables, 27 constraints
        
        int order = 3;//what order of approximation to use
        
        DoubleMatrix1D derivVec = getCoordDerivs(block);
        
        RotMatrix r = new RotMatrix();
        
        double ans[][] = new double[9][];
        for(int at=0; at<9; at++)
            ans[at] = r.add( blockActualCoordOld[block][at], derivVec.viewPart(3*at,3).toArray() );
        
        DoubleMatrix1D secondDerivVec=null, thirdDerivVec=null;
        
        if(order>1){
            secondDerivVec = getCoordSecondDerivs(block,derivVec);
            secondDerivVec.assign( Functions.div(2.) );//need 1/2! factor for Taylor series
            for(int at=0; at<9; at++)
                ans[at] = r.add( ans[at], secondDerivVec.viewPart(3*at,3).toArray() );
            
            if(order>2){
                thirdDerivVec = getCoordThirdDerivs(block,derivVec,secondDerivVec);
                thirdDerivVec.assign( Functions.div(6.) );//need 1/2! factor for Taylor series
                for(int at=0; at<9; at++)
                    ans[at] = r.add( ans[at], thirdDerivVec.viewPart(3*at,3).toArray() );
            }
        }
        
        
        //DEBUG!!!
        //track steps in progression from t=0 (old) to t=1 (new).  Should start with right derivatives
        /*
        double tSteps[] = new double[] {0,0.001,0.002,0.003,0.005,0.008,0.01,0.02,0.03,0.04,0.05,0.08,0.1,0.2,0.5
                ,0.8,1};
        
        System.out.println("t firstCA_x lastCA_x CAC1 CAC2 NC2 CACA1");
        
        for(double t : tSteps){
            
            double coords[][] = new double[9][];
            for(int at=0; at<9; at++)
                coords[at] = r.add( blockActualCoordOld[block][at], derivVec.viewPart(3*at,3).copy().
                        assign(Functions.mult(t)).toArray() );
            
            if(order>1){
                for(int at=0; at<9; at++)
                    coords[at] = r.add( coords[at], secondDerivVec.viewPart(3*at,3).copy().
                            assign(Functions.mult(t*t/2)).toArray() );
                
                if(order>2){
                    for(int at=0; at<9; at++)
                        coords[at] = r.add( coords[at], thirdDerivVec.viewPart(3*at,3).copy().
                            assign(Functions.mult(t*t*t/6)).toArray() );
                }
            }
            
            //some dists
            double cac1 = r.norm(r.subtract(coords[1],coords[2]));
            double cac2 = r.norm(r.subtract(coords[4],coords[5]));
            double nc2 = r.norm(r.subtract(coords[3],coords[5]));
            double caca1 = r.norm(r.subtract(coords[1],coords[4]));

            System.out.println(t+" "+coords[1][0]+" "+coords[7][0]+" "+cac1+" "+cac2+" "+
                    nc2+" "+caca1);
        }
        
        
        int qqq = 0;
        */
        
        
        
        
        return ans;
    }
    
    //out approximation is based on the anchor atoms moving at constant velocity
    //from the blockActualCoordOld to the blockActualCoordNew positions
    DoubleMatrix1D getCoordDerivs(int block){
        
        DoubleMatrix2D M = DoubleFactory2D.dense.make(27,27);
        DoubleMatrix2D C = DoubleFactory2D.dense.make(27,1);//right-hand side of lin eq.
        
        
        int count = 0;//constraint counter

        //for convenience, free-atom coords (first and last 2 of blockAtoms) added as lin constraints
        //CAN REARRANGE FOR SPEEDUP IF NEEDED
        for(int freeAt : new int[] {0,1,7,8}){
            for(int dim=0; dim<3; dim++){
                M.set(count,3*freeAt+dim,1);
                C.set(count,0,blockActualCoordNew[block][freeAt][dim]-blockActualCoordOld[block][freeAt][dim]);
                count++;
            }
        }
        
        //distance constraints
        for(int at=1; at<=6; at++){//all 1,2-distances involving atoms with undetermined positions
            makeDistConstr(at,at+1,count,M,block);
            count++;
        }
        
        
        //1,3-distance constraints (proxy for angle constr)
        for(int at=0; at<=6; at++){//all 1,3-distances involving atoms with undetermined positions
            makeDistConstr(at,at+2,count,M,block);
            count++;
        }
        
        
        //CA-CA distance constraints (proxy for omega constr)
        for(int at : new int[] {1,4}){//These are the first and second CAs
            makeDistConstr(at,at+3,count,M,block);
            count++;
        }
        
        
        //now solve the linear system
        return Algebra.DEFAULT.solve(M,C).viewColumn(0);
    }
    
    
    
    
    
    void makeDistConstr(int at1, int at2, int row, DoubleMatrix2D M, int block){
        //we're building a linear system Mx = C where x are the time derivatives of the atoms
        //in the given closure block.  
        //We want to constrain the time derivative of the distance between atoms at1 and at2
        //(numbered within block) to be 0, and put this constraint into the given row of M
        //right-hand side of linear system stays 0
        
        //ok if the atoms' coords are x and y
        //then our constraint is just (dx/dt - dy/dt) dot (x-y) = 0
        //we'll use the pre-motion x and y for this (since we know them:
        //blockActualCoordNew is not physically meaningful for atoms inside block
        //but just puts the atoms in position to be put in place by an LCA)
        for(int dim=0; dim<3; dim++){
            double diff = blockActualCoordOld[block][at1][dim] - blockActualCoordOld[block][at2][dim];
            M.set(row, 3*at1+dim, diff);
            M.set(row, 3*at2+dim, -diff);
        } 
    }
    
    
    
    DoubleMatrix1D getCoordSecondDerivs(int block, DoubleMatrix1D coordDerivs){
        //same idea as getCoordDerivs, but getting second derivs
        
        DoubleMatrix2D M = DoubleFactory2D.dense.make(27,27);
        DoubleMatrix2D C = DoubleFactory2D.dense.make(27,1);//right-hand side of lin eq.
        
        
        int count = 0;//constraint counter

        //for convenience, free-atom coords (first and last 2 of blockAtoms) added as lin constraints
        //CAN REARRANGE FOR SPEEDUP IF NEEDED
        for(int freeAt : new int[] {0,1,7,8}){
            for(int dim=0; dim<3; dim++){
                M.set(count,3*freeAt+dim,1);
                //second derivatives of free-atom coords held to 0: constant-velocity motion
                count++;
            }
        }
        
        //distance constraints
        for(int at=1; at<=6; at++){//all 1,2-distances involving atoms with undetermined positions
            makeDistConstr2(at,at+1,count,M,C,block,coordDerivs);
            count++;
        }
        
        
        //1,3-distance constraints (proxy for angle constr)
        for(int at=0; at<=6; at++){//all 1,3-distances involving atoms with undetermined positions
            makeDistConstr2(at,at+2,count,M,C,block,coordDerivs);
            count++;
        }
        
        
        //CA-CA distance constraints (proxy for omega constr)
        for(int at : new int[] {1,4}){//These are the first and second CAs
            makeDistConstr2(at,at+3,count,M,C,block,coordDerivs);
            count++;
        }
        
        
        //now solve the linear system
        return Algebra.DEFAULT.solve(M,C).viewColumn(0);
    }
    
    
    void makeDistConstr2(int at1, int at2, int row, DoubleMatrix2D M, DoubleMatrix2D C, int block,
            DoubleMatrix1D derivs1){
        //we're building a linear system Mx = C where x are the 2nd time derivatives of the atoms
        //in the given closure block.  
        //We want to constrain the 2nd time derivative of the distance between atoms at1 and at2
        //(numbered within block) to be 0, and put this constraint into the given row of M
        //right-hand side of linear system stays 0
        
        //ok if the atoms' coords are x and y
        //constraint is (x-y) dot d^2(x-y)/dt^2 = - d(x-y)/dt dot d(x-y)/ dt
        //we'll use the pre-motion x and y for this (since we know them:
        //blockActualCoordNew is not physically meaningful for atoms inside block
        //but just puts the atoms in position to be put in place by an LCA)
        for(int dim=0; dim<3; dim++){
            double diff = blockActualCoordOld[block][at1][dim] - blockActualCoordOld[block][at2][dim];
            M.set(row, 3*at1+dim, diff);
            M.set(row, 3*at2+dim, -diff);
        }
        
        DoubleMatrix1D diffDeriv1 = derivs1.viewPart(3*at1,3).copy().assign( 
                derivs1.viewPart(3*at2,3), Functions.minus );
        
        C.set(row, 0, -diffDeriv1.zDotProduct(diffDeriv1) );
    }
    
    
    
    
    DoubleMatrix1D getCoordThirdDerivs(int block, DoubleMatrix1D coordDerivs, DoubleMatrix1D coordDerivs2){
        //same idea as getCoordDerivs, but getting second derivs
        
        DoubleMatrix2D M = DoubleFactory2D.dense.make(27,27);
        DoubleMatrix2D C = DoubleFactory2D.dense.make(27,1);//right-hand side of lin eq.
        
        
        int count = 0;//constraint counter

        //for convenience, free-atom coords (first and last 2 of blockAtoms) added as lin constraints
        //CAN REARRANGE FOR SPEEDUP IF NEEDED
        for(int freeAt : new int[] {0,1,7,8}){
            for(int dim=0; dim<3; dim++){
                M.set(count,3*freeAt+dim,1);
                //second derivatives of free-atom coords held to 0: constant-velocity motion
                count++;
            }
        }
        
        //distance constraints
        for(int at=1; at<=6; at++){//all 1,2-distances involving atoms with undetermined positions
            makeDistConstr3(at,at+1,count,M,C,block,coordDerivs,coordDerivs2);
            count++;
        }
        
        
        //1,3-distance constraints (proxy for angle constr)
        for(int at=0; at<=6; at++){//all 1,3-distances involving atoms with undetermined positions
            makeDistConstr3(at,at+2,count,M,C,block,coordDerivs,coordDerivs2);
            count++;
        }
        
        
        //CA-CA distance constraints (proxy for omega constr)
        for(int at : new int[] {1,4}){//These are the first and second CAs
            makeDistConstr3(at,at+3,count,M,C,block,coordDerivs,coordDerivs2);
            count++;
        }
        
        
        //now solve the linear system
        return Algebra.DEFAULT.solve(M,C).viewColumn(0);
    }
    
    
    void makeDistConstr3(int at1, int at2, int row, DoubleMatrix2D M, DoubleMatrix2D C, int block,
            DoubleMatrix1D derivs1, DoubleMatrix1D derivs2){
        //third derivatives
        //constraint is (x-y) dot d^3(x-y)/dt^3 = - 3 d(x-y)/dt dot d^2(x-y)/ dt^2
        for(int dim=0; dim<3; dim++){
            double diff = blockActualCoordOld[block][at1][dim] - blockActualCoordOld[block][at2][dim];
            M.set(row, 3*at1+dim, diff);
            M.set(row, 3*at2+dim, -diff);
        }
        
        DoubleMatrix1D diffDeriv1 = derivs1.viewPart(3*at1,3).copy().assign( 
                derivs1.viewPart(3*at2,3), Functions.minus );
        
        DoubleMatrix1D diffDeriv2 = derivs2.viewPart(3*at1,3).copy().assign( 
                derivs2.viewPart(3*at2,3), Functions.minus );
        
        C.set(row, 0, -3*diffDeriv1.zDotProduct(diffDeriv2) );
    }
    
    
    
    double getBlockLinApproxDeviation(int block, int soln, double[][] linApprox){
        //sum of square deviations of linear approximation from closure block
        //from the given solution in blockLCAs[block]
        //return inf if solution violates maxCATrans constraints
        
        //THIS CAN ALSO BE OPTIMIZED IF NEEDED
        //BUT THE BEST WAY TO HANDLE OPTIMIZATION MIGHT BE TO ITERATIVELY GET CLOSURE
        //BY MOVING CA FROM LIN APPROX ANYWAY
        
        double[][][] rotMatrices = blockLCAs[block].rotMatrices[0][soln];
        double ssd = 0;
        RotMatrix r = new RotMatrix();
        
        //deviation limited to the five middle blockAtoms (#2-6)
        
        //#2-4 (up to middle CA) suffer the first rotation, about first CA (atom #1)
        for(int at=2; at<=4; at++){
            double[] solnCoord = r.rotateAboutPt(blockActualCoordNew[block][at], rotMatrices[0], 
                    blockActualCoordNew[block][1]);//atom coords based on LCA soln
            
            if(at==4){//maxCATrans filter
                double origCoord[] = m.residue[closureBlock[block]].getAtomByName("CA").coord;
                for(int dim=0; dim<3; dim++){
                    if( Math.abs(origCoord[dim]-solnCoord[dim]) > maxCATrans )
                        return Double.POSITIVE_INFINITY;
                }
            }
            
            ssd += r.normsq( r.subtract(solnCoord,linApprox[at]) );
        }
        
        
        //#5-6 (middle C', last N) suffer second rotation, about last CA (atom #7)
        for(int at=5; at<=6; at++){
            double[] solnCoord = r.rotateAboutPt(blockActualCoordNew[block][at], rotMatrices[1], 
                    blockActualCoordNew[block][7]);//atom coords based on LCA soln
            
            ssd += r.normsq( r.subtract(solnCoord,linApprox[at]) );
        }
        
        return ssd;
    }
    
    
    
    
    void recordBlockActualCoords(boolean useNew){
        //For each closure block, record actual coordinates for the mainchain atoms involved
        
        //decide which field to record in
        double[][][] curField = blockActualCoordOld;
        if(useNew)
            curField = blockActualCoordNew;
        
        
        for(int block=0; block<closureBlock.length; block++){
            //we record the 9 mainchain atoms for the 3 residues involved in the block
            //(closureBlock[block]-1 through closureBlock[block+1])
            //record atoms in chain order
            
            
            for(int res=0; res<3; res++){
                Residue curRes = m.residue[closureBlock[block]+res-1];
                curField[block][3*res] = m.getActualCoord(curRes.getAtomNameToMolnum("N"));
                curField[block][3*res+1] = m.getActualCoord(curRes.getAtomNameToMolnum("CA"));
                curField[block][3*res+2] = m.getActualCoord(curRes.getAtomNameToMolnum("C"));
            }
        }
    }
    
    
    /*

    void calcSolns(int blockNum){
        //Based on the anchor atom coords for the specified closure block currently specified
        //in the actualCoordinates, get the solutions (for the middle CA rot parameter tan(phi/2) )
        //and put them in blockSolns

        DoubleMatrix1D N0, CA2, C2;//make these relative to CA1, with CA2 on z-axis and CA1
        //starting on x axis

        MultivariatePoly CA1Polys[] = new MultivariatePoly[3];//express Ca1
        CA1Polys[0] = a*cosphi;
        //then asinphi and b
        
        
        //generate the linear functions of (cos phi,sin phi) that describe the N1 and C1 coords
        //in our local coord. system of course
        //they are expressed as MultivariatePolys
        //so coord i of N1 is (N1Polys[i]/N1Polys[6] + N1Polys[i+3]).  Same for C1.
        MultivariatePoly N1Polys[] = new MultivariatePoly[7];
        MultivariatePoly C1Polys[] = new MultivariatePoly[7];
        
        getPepPlanePolys(CA1Polys,0,N0,N1Polys);
        getPepPlanePolys(CA1Polys,CA2,C2,C1Polys);
        
        MultivariatePoly N1denom = N1Polys[6];
        MultivariatePoly C1denom = C1Polys[6];
        
        //Let's express our equation in the form Q1/(N1denom*C1denom) + Q2/N1denom + Q3/C1denom + Q4 = 0
        //where the Q's are quadratic in cos phi, sin phi
        MultivariatePoly Q1=0, Q2=0, Q3=0, Q4=0;
        
        //our equation is (N1-CA1) dot (C1-CA1) = (what it should be for a residue);

        for(int dim=0; dim<3; dim++){//dimensions of the dot product
            Q1.add(N1Polys[dim].times(C1Polys[dim]));
            
            Q2.add(N1Polys[dim].times(C1Polys[dim+3]));
            Q2.subtract(N1Polys[dim].times(CA1Polys[dim]));
            
            Q3.add(C1Polys[dim].times(C1Polys[dim+3]));
            Q3.subtract(C1Polys[dim].times(Ca1Polys[dim]));
            
            Q4.add(N1Polys[dim+3].times(C1Polys[dim+3]));
            Q4.subtract(N1Polys[dim+3].times(CA1Polys[dim]));
            Q4.subtract(C1Polys[dim+3].times(CA1Polys[dim]));
            Q4.add(CA1Polys[dim].times(CA1Polys[dim]));
            
            Q4 -= getTargetDotProduct("N",closureBlocks[blockNum],"CA",closureBlocks[blockNum],"C",closureBlocks[blockNum]);
        }
        
                
        
        //Now we can build our quartic (total degree) w.r.t. cos phi and sin phi
        //we must multiply through by (N1denom*C1denom)
        MultivariatePoly quarticPoly = Q1.copy();
        quarticPoly.add(Q2.times(C1denom));
        quarticPoly.add(Q3.times(N1denom));
        quarticPoly.add(Q4.times(C1denom).times(N1denom));
        
        
        //Finally, convert to a univariate octic polynomial in tan(phi/2)
        //and solve
        double[] octic = quarticPoly.halfAngTangentForm();
        
        //ideally we'll limit solutions to our desired range...
        blockSolns[blockNum] = SturmSolver.solveOctic(octic);
    }

    
    
    void getPepPlanePolys(MultivariatePoly[] CA1Polys, DoubleMatrix1D anchor, DoubleMatrix1D beyondAnchor,
            MultivariatePoly[] ans){
        //Let CA1Polys express the middle-CA coords of a closure block
        //we have (local-system) coords of the anchor CA and the atom beyond it
        //express the atom closest to that anchor of the middle residue in terms of polynomials ans: 
        //so coord i of N1 is (N1Polys[i]/N1Polys[6] + N1Polys[i+3]).  Same for C1.
        
        //First, we want to express the unknown atom closest to the anchor as a vector z=[x,m1*x+b1,m2*x+b2]
        //this atom forms a triangle of known shape with anchor and beyondAnchor
        //we get two linear equations based on the dot products of edges of this triangle
        double dp1 = getTargetDotProduct(...);//(z-anchor)dot(beyondAnchor-anchor)=dp1
        
        //our equations are u dot z = g, v dot z = h
        DoubleMatrix1D u = b;
        
        
        
        
    }
    
    
    double getTargetDotProduct(String name1, int res1, String name2, int res2, String name3, int res3){
        //When trying to close a closure block, this is our target for (atom 1-atom2) dot (atom3-atom2)
        //Each atom specified as name, molResNum
        //We'll currently do this based on the original (Atom.coord) coordinates for the atoms
        //with the goal of being able to exactly reproduce the starting geometry
        //but in some cases one might instead want the dot product based on ideal bond lengths, angles, and omega
        //and this function should only be used for cases where the dot product could be obtained based 
        //on those ideal values (i.e. the dot product is not changed by any "soft" DOFs)
        double[] at1 = m.residue[res1].getAtomByName(name1).coord;
        double[] at2 = m.residue[res2].getAtomByName(name2).coord;
        double[] at3 = m.residue[res3].getAtomByName(name3).coord;
        
        double ans = 0;
        for(int dim=0; dim<3; dim++)
            ans += (at1[dim]-at2[dim])*(at3[dim]-at2[dim]);
        
        return ans;
    }
    */
    
    
    
    void applySoln(int blockNum, int solnNum){
        
        if(solnNum<-1 || solnNum>=blockNumSolns[blockNum])
            throw new RuntimeException("ERROR: Closure block "+blockNum+" has no solution # "+solnNum);
        
        blockLCAs[blockNum].doPerturbationMotion(solnNum);
    }
    
    
    void tryAltClosures(EnergyFunction efunc){
        //try alternate closures for each closure block at a time
        //picking each to minimize efunc
        
        for(int res=startRes; res<=endRes; res++)
            genChi1[res-startRes] = Perturbation.getGenChi1(m, res);
        
        
        for(int block=0; block<closureBlock.length; block++){
            
            calcSolns(block);
            
            int bestSoln = -1;
            double bestE = Double.POSITIVE_INFINITY;

            double curCoords[][] = new double[3][];
            for(int a=0; a<3; a++)
                curCoords[a] = storeResActualCoord(closureBlock[block]+a-1);
            
            for(int soln=0; soln<blockNumSolns[block]; soln++){
              
                //first check if middle CA moved too much
                boolean movedTooMuch = false;
                RotMatrix r = new RotMatrix();
                
                double firstCA[] = m.getActualCoord( m.residue[closureBlock[block]-1].getAtomNameToMolnum("CA") );
                double middleCA[] = m.getActualCoord( m.residue[closureBlock[block]].getAtomNameToMolnum("CA") );//middle CA before rotation
                double middleCAOrig[] = m.residue[closureBlock[block]].getAtomByName("CA").coord;
                double middleCASoln[] = r.rotateAboutPt(middleCA, blockLCAs[block].rotMatrices[0][soln][0],
                    firstCA);//middle CA coords based on LCA soln
                
                for(int dim=0; dim<3; dim++){
                    if( Math.abs(middleCAOrig[dim]-middleCASoln[dim]) > maxCATrans )
                        movedTooMuch = true;
                }

                
                if(!movedTooMuch){
                    
                    applySoln(block, soln);
                    double E = Double.POSITIVE_INFINITY;
                    boolean SCOK = true;//can get good sidechains (i.e. can close any proline rings)
                    
                    for(int res=closureBlock[block]-1; res<=closureBlock[block]+1; res++)
                        SCOK = SCOK && fixSC(res);
                    
                    if(SCOK)
                        E = efunc.getEnergy(0);
                    
                    if(E<bestE){
                        bestE = E;
                        bestSoln = soln;
                    }
                    
                    for(int a=0; a<3; a++)
                        restoreResActualCoord(closureBlock[block]+a-1,curCoords[a]);
                }
            }
            
            if(bestSoln==-1)
                System.out.println("Warning: Trying alternate closures but can't even get the original loop closure.");
            
            
            applySoln(block,bestSoln);
            
            for(int res=closureBlock[block]-1; res<=closureBlock[block]+1; res++)
                fixSC(res);//we know this will work because we already tried it in the loop above
        }
        
    }
    
    
    boolean fixSC(int res){
        //orient the sidechain properly, for the given residue (numbered within this FBFB)
        //if we can't (this should only happen for proline, then return false and we'll revert the backbone
        boolean outcome = m.idealizeResSidechain(m.residue[res]);
        if(!outcome)
            return false;
        
        Perturbation.setGenChi1(m, res, genChi1[res-startRes]);
        return true;
    }
    
    
    
    //In CCDMinimizer, we can first setFreeDOFs
            //and then try out phi solns one at a time, to get to the min (discrete version of coord descent)
       
    
    
    //support the gradient, using alternating peptide planes to parameterize
    //and thus to identify our solutions of interest?
    //BETTER TO USE LINEAR APPROX RATHER THAN ALTERNATING PEP PLANES

    

    //information on parameters
    double getMeshWidth(int paramNum){
        //These values are meant to be similar to those in DegreeOfFreedom.getMeshWidth
        if(freeDOFs[paramNum].type==ANCHORDIH)
            return 1.;
        else if(freeDOFs[paramNum].type==FIRSTCAADJ)
            return 0.1;
        else//peptide plane rotations
            return 1.;
    }
    
    boolean isParamAngle(int paramNum){
        if(freeDOFs[paramNum].type==FIRSTCAADJ)
            return false;
        else//either dihedral of peptide plane rotation angle
            return true;
    }
    
    int[] paramAffectedRes(int paramNum){
        //molecule residue numbers for all residues
        //whose conformations are dependent on the value of the given parameter
        
        BBFreeDOF bfd = freeDOFs[paramNum];
        
        if(bfd.type==ANCHORDIH){
            if(bfd.anchorDihNum<=1)//everything from s†ar†Res †hrough firs† closure block
                return new int[] {startRes,startRes+1,startRes+2,startRes+3};
            else
                return new int[] {endRes-3,endRes-2,endRes-1,endRes};
        }
        else {//free peptide plane parameters: include the free peptide plane
            //and the closure blocks on either side
            int plane = freePepPlanes[bfd.freePlaneNum];
            return new int[] {plane-2,plane-1,plane,plane+1,plane+2,plane+3};
        }
    }
            

}