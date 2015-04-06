
import cern.colt.matrix.DoubleMatrix1D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
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
//	QuantumChangeEFcn.java
//
//	Version:           2.2 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class QuantumChangeEFcn extends EnergyFunction {
    //this is an initial test of EPIC's ability to represent quantum energies
    //similar idea to PBChangeEFcn, except we will be applying the dihedrals to a separate molecule
    //hard-coding things to aspartame bc its just a test
     
    
    
    DoubleMatrix1D DOFVals;
    DoubleMatrix1D origDOFVals;
    
    int mainSystem;//system for energy change: either 0=ASP, 1=PHE 41, or 2=both (pairwise)
    
    double origE;//value at starting coordinates
    
    Molecule aspartame;//the molecule we'll be messing with
    
    
    static String NWChemFolder="OSPREY_nwchem";
    
    EnergyFunction oldEF = null;//the energy function this is replacing
    //used to make SVE
    
    
    
    public QuantumChangeEFcn(DoubleMatrix1D x, int res1, int res2, EnergyFunction oldEF){
        
        //load aspartame w/o autofixing
        boolean au = EnvironmentVars.autoFix;
        EnvironmentVars.autoFix = false;
        try {
            FileInputStream is = new FileInputStream(NWChemFolder+"/aspartame.pdb");
            aspartame = new Molecule();
            new PDBChemModel(aspartame, is);
        }
        catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
        EnvironmentVars.autoFix = au;
        
        
        
        origDOFVals = x.copy();
        DOFVals = x;//changes will be made to this elsewhere, which will be read off here
        
        if(res2==-1)
            mainSystem = res1;
        else//pairwise
            mainSystem = 2;
        
        //set origE
        applyDOFVals(origDOFVals);
        origE = getSystemE();
        
        this.oldEF = oldEF;
    }
    
    
    public double getEnergy(){
        
        if(mainSystem==2){
            
            applyDOFVals(DOFVals);
            double E = getSystemE();
            
            applyDOFVal(2,origDOFVals.get(2));//use original Phe dihedrals
            applyDOFVal(3,origDOFVals.get(3));
            E -= getSystemE();
            
            applyDOFVal(0,origDOFVals.get(0));//use original Asp dihedrals
            applyDOFVal(1,origDOFVals.get(1));
            applyDOFVal(2,DOFVals.get(2));
            applyDOFVal(3,DOFVals.get(3));
            E -= getSystemE();
            
            E += origE;
            return E;
        }
        else{
            applyDOFVals(DOFVals);
            return getSystemE()-origE;
        }
    }
   
    
    
    void applyDOFVals(DoubleMatrix1D x){
        //apply DOF values for the main system
        if(mainSystem==0 || mainSystem==2){
            applyDOFVal(0,x.get(0));
            applyDOFVal(1,x.get(1));
        }
        if(mainSystem==1){
            applyDOFVal(2,x.get(0));
            applyDOFVal(3,x.get(1));
        }
        if(mainSystem==2){
            applyDOFVal(2,x.get(2));
            applyDOFVal(3,x.get(3));
        }
    }
    
    
    void applyDOFVal(int dihNum, double chi){
        //hard-coding the 4 sidechain dihedrals for aspartame, indexed by dihNum...
        if(dihNum==0)
            aspartame.setTorsion(0,1,2,3,chi,new int[] {4,5,11,12},4,false);//asp chi1
        if(dihNum==1)
            aspartame.setTorsion(1,2,3,4,chi,new int[] {5},1,false);//asp chi2
        if(dihNum==2)
            aspartame.setTorsion(14,15,16,17,chi,new int[] {18,19,20,21,22,26,27,29,30,31,32,33},12,false);//phe chi1
        if(dihNum==3)
            aspartame.setTorsion(15,16,17,18,chi,new int[] {19,20,21,22,29,30,31,32,33},9,false);//phe chi2
    }
    
    
    double getSystemE(){
        //get energy of aspartame in its current conformation
        double E = 0;
        //new KSParser().saveMolecule(aspartame, NWChemFolder+"/molec.pdb", 0);
        //we need more precision...the little bond length/angle fluctuations from PDB roundoff 
        //seem to dominate the dihedral-related energy changes
        
        for(int tryNum=0; true; tryNum++){
            
            try{
                //writeNWFile();
                BufferedWriter bw=new BufferedWriter(new FileWriter(NWChemFolder+"/molec.coords.txt"));
                for(int a=0; a<aspartame.numberOfAtoms; a++){
                    double x[] = aspartame.getActualCoord(a);
                    bw.append( "\t"+aspartame.atom[a].elementType+
                            " "+x[0]+" "+x[1]+" "+x[2] );
                    bw.newLine();
                }

                bw.close();


                Process p = Runtime.getRuntime().exec(NWChemFolder+"/getE_SCF");
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                E = Double.valueOf( br.readLine().trim() );
                E *= 627.509469;//convert from hartrees (NWChem unit) to kcal/mol (OSPREY unit)
                br.close();
                break;
            }
            catch(Exception e){
                
                if(tryNum>5)
                    throw new RuntimeException("ERROR evaluating quantum energy, even with retrying "+e.getMessage());
                else {
                    System.out.println("Warning: Needing to retry quantum energy: "+e.getMessage());
                }
            }
        }

        return E;
    }
    
    
    
    void writeNWFile() throws Exception {
        BufferedWriter bw=new BufferedWriter(new FileWriter(NWChemFolder+"/coords.txt"));
        for(int i=0; i<aspartame.numberOfAtoms; i++){
            bw.append("\t"+aspartame.atom[i].elementType+" "+aspartame.actualCoordinates[3*i]+" "+
                    aspartame.actualCoordinates[3*i+1]+" "+aspartame.actualCoordinates[3*i+2]);
        }
        bw.close();
    }

    
    
    
   
    @Override
    public void setupPartialComputation(int[][] residues) {
        //not doing partial comp...would require PB perturbation theory probably
    }

    @Override
    public double getEnergy(int a) {
        return getEnergy();
    }
    
    
    
    /*@Override
    public Amber96ext getAmber96ext(){
        if(ae!=null)
            return ae;
        throw new RuntimeException("ERROR: no Amber96ext assigned to QuantumChangeEFcn");
    }*/
    
}
