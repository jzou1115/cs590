
import java.io.BufferedReader;
import java.io.InputStreamReader;
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
//	PBChangeEFcn.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class PBChangeEFcn extends EnergyFunction {
    //This energy function is meant to capture continuous variation in Poisson-Boltzmann energies due to dihedral angle changes
    //(so the molecule involved is not expected to have mutations)
    //we start with a structure (ideal rotamers)
    //then for intra+shell energies, we return Poisson-Boltzmann energy of the molecule w/ changed dihedrals (expected to be at one res) - orig E
    //for pairwise, we do total PB E - (E with only res1 changed) - (E with only res2 changed) + orig E
    //since this is just a first investigation into representing Poisson-Boltzmann energies,
    //it is not particularly speed-optimized.  PB calc will likely be the bottleneck anyway.  
    
    Molecule m;//The molecule whose conformational changes we will be calculating energies for
    Molecule mCopy;//A copy of the molecule used for outputting PDB files, which Delphi will read
    double origE;
    
    boolean pairwise = false;//is this a pairwise energy?
    //if so here are the molecule residue numbers for the pair:
    int res1;
    int res2;
    
    static String delphiFolder = "OSPREY_delphi";
    
    public PBChangeEFcn(Molecule molec, int str1, int strResNum1, int str2, int strResNum2){
        //if pairwise, then we must have valid strand and strand res. num. arguments
        //non-pairwise flagged by str2==-1 (e.g. this is found for shell runs in RotamerSearch)
        
        m = molec;
        
        try{
            mCopy = (Molecule)KSParser.deepCopy(molec);
        }
        catch(Exception e){
            throw new RuntimeException("ERROR deepCopying molec in new PBChangeEFcn()");
        }
        
        mCopy.backupCoordinates = new double[mCopy.actualCoordinates.length];
        fullCopy(mCopy.actualCoordinates,mCopy.backupCoordinates);
        
        origE = mCopyE();
        
        if(str2!=-1){//pairwise
            pairwise = true;
            res1 = m.strand[str1].residue[strResNum1].moleculeResidueNumber;
            res2 = m.strand[str2].residue[strResNum2].moleculeResidueNumber;
        }
    }

    
    @Override
    public double getEnergy() {
        if(pairwise){
            fullCopy(m.actualCoordinates,mCopy.actualCoordinates);
            double ans = mCopyE();
            ans -= singleResE(res1);
            ans -= singleResE(res2);
            ans += origE;
            return ans;
        }
        else {
            //copy over coordinate and evaluate energy
            fullCopy(m.actualCoordinates,mCopy.actualCoordinates);
            return mCopyE()-origE;
        }
    }
    
    double singleResE(int molResNum){
        //create a hybrid molecule that has coordinates of m at molResNum, original (backed-up) coordinates elsewhere
        fullCopy(mCopy.backupCoordinates,mCopy.actualCoordinates);
        Residue res = m.residue[molResNum];
        for(Atom at : res.atom){
            int offset = 3*at.moleculeAtomNumber;
            System.arraycopy(m.actualCoordinates, offset, mCopy.actualCoordinates, offset, 3);
        }
        return mCopyE();   
    }
    
    
    void fullCopy(double[] a, double[] b){
        //full copy of double array
        System.arraycopy(a,0,b,0,a.length);
    }
    
    
    double mCopyE(){
        //get the current PB energy of mCopy
        //mCopy doesn't use Atom.coord except when saving the molecule, so we can just output actualCoordinates using saveMolec
        //without messing anything up
        
        //save the molecule to a PDB file and run Delphi on it
        new KSParser().saveMolecule(mCopy, delphiFolder+"/struct.pdb", 0);
        double E = 0;
        
        try{
            Process p = Runtime.getRuntime().exec(delphiFolder+"/getE");
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            E = Double.valueOf( br.readLine().trim() );
            E *= RotamerSearch.constRT;//convert from thermal energies (Delphi unit) to kcal/mol (OSPREY unit)
            br.close();
        }
        catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        return E;
    }
    
    
    

    @Override
    public void setupPartialComputation(int[][] residues) {
        //not doing partial comp...would require PB perturbation theory probably
    }

    @Override
    public double getEnergy(int a) {
        return getEnergy();
    }
    
    
    
}
