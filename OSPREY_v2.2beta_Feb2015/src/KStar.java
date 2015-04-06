
import java.io.IOException;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;

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

////////////////////////////////////////////////////////////////////////////////////////////
// KStar.java
//
//  Version:           2.2 beta
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
//    PGC        Pablo Gainza C.       Duke University         pablo.gainza@duke.edu
////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2009)
 *	
 */


/**
 * This is the main class for the KStar program; essentially just a wrapper for the KSParser class.
 */
public class KStar
{
    	
	public static void main (String[] args) throws IOException
	{   
            
            /*
             //PG9 graft
            EnvironmentVars.autoFix = false;//so we don't need to get templates
            
            RamachandranChecker.getInstance().readInputFiles( new String[] { 
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-gly-sym.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-pro.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-general.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-prepro.data" } );
            
            //Create StructGrafter for structures struct1, struct2
            StructGrafter sg = new StructGrafter("../input_structs/4LSS.renum.pdb",
                    "../input_structs/3U2S.renum.pdb");
            
            //align struct2 to struct1
            sg.SSAlignment(new String[][] {new String[]{"791","797"},new String[]{"805","812"}}, 
                    new String[][] {new String[]{"991","997"},new String[]{"1005","1012"}},
                    "4LSS.3U2S.ssalign.pdb");
            //OR anchorAlignment
            //sg.anchorAlignment("797","802","97","102","4LSS.4JM2.anchoralign.pdb");
            
            //graft structure together
            //start with whole antibody part of 4LSS
            sg.initCompositeStructure("701","1616");
            
            //graft in loop
            sg.graftLoop("795", "803", "995", "1003");//was 97-105...
            //add ligand
            sg.transferLigand("120","573");
            
            //let's see what the structure is like now
            sg.struct3.connectivity12Valid = true;//let's play pretend first
            new KSParser().saveMolecule(sg.struct3, "4LSS.3U2S.open.pdb", 0);
            //These settings give a good 4LSS/4JM2 graft geometry- and BB clash-wise.  
            
            //close loops
            //replaced 794-804 (PG9_graft)
            //PG9_graft2 did 796-802
            //PG9_graft3 will do 795-803
            //sending graft to mhall44, graft2 to mah43, and graft3 to mark.hallen
            sg.matchLoopAnchors("795","800","800L","803","4LSS.3U2S.grafted.pdb");
            
            
            
            //make graft with CD4bs ligand ("graftc")
            StructGrafter sg2 = new StructGrafter("4LSS.3U2S.grafted.pdb","../input_structs/4LSS.renum.pdb");
            sg2.initCompositeStructure("701","1616");
            sg2.transferLigand("44","492");
            sg2.struct3.connectivity12Valid = true;
            new KSParser().saveMolecule(sg2.struct3, "4LSS.3U2S.graftc.pdb", 0);

            System.exit(0);
            */
            
            /*
            //CH103 epitope graft
            EnvironmentVars.autoFix = false;//so we don't need to get templates
            
            RamachandranChecker.getInstance().readInputFiles( new String[] { 
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-gly-sym.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-pro.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-general.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-prepro.data" } );
            
            //Create StructGrafter for structures struct1, struct2
            StructGrafter sg = new StructGrafter("../input_structs/4LSS.renum.pdb",
                    "../input_structs/4JAN.renum.pdb");
            
            //align struct2 to struct1
            sg.SSAlignment(new String[][] {new String[]{"797","797"},new String[]{"801","801"}}, 
                    new String[][] {new String[]{"797","797"},new String[]{"801","801"}},
                    "4LSS.4JAN.ssalign.pdb");
            
            //tried aligning 797-7 and 801-12; this gave clashes around 800a, because.
            //tried aligning just 797 and 801 but this causes clashes elsewhere, e.g. gp120-light chain
            //seems there's a substantial twist between the VRC01 and CH103 CDRH3 loops
            //(in the secondary-structure frame of reference)
            //and, because there's a fairly large area of contact between CH103 and gp120,
            //tilting the gp120 to match this twist (as is needed to accommodate the loop without local clashes)
            //causes gp120 to clash elsewhere
            
            
            //graft structure together
            //start with whole antibody part of 4LSS
            sg.initCompositeStructure("701","1616");
            
            //graft in loop
            sg.graftLoop("794", "802", "794", "802");//was 97-105...
            //add ligand
            sg.transferLigand("255","605");
            
            //let's see what the structure is like now
            sg.struct3.connectivity12Valid = true;//let's play pretend first
            new KSParser().saveMolecule(sg.struct3, "4LSS.4JAN.open.pdb", 0);
            //These settings give a good 4LSS/4JM2 graft geometry- and BB clash-wise.  
            
            //close loops
            sg.matchLoopAnchors("794","800","800B","802","4LSS.4JAN.grafted.pdb");
            
            
            
            //make graft with CD4bs ligand ("graftc")
            StructGrafter sg2 = new StructGrafter("4LSS.4JAN.grafted.pdb","../input_structs/4LSS.renum.pdb");
            sg2.initCompositeStructure("701","1616");
            sg2.transferLigand("44","492");
            sg2.struct3.connectivity12Valid = true;
            new KSParser().saveMolecule(sg2.struct3, "4LSS.4JAN.graftc.pdb", 0);

            System.exit(0);
            */
            
            /*
             //PGT135 graft
            EnvironmentVars.autoFix = false;//so we don't need to get templates
            
            RamachandranChecker.getInstance().readInputFiles( new String[] { 
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-gly-sym.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-pro.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-general.data",
                "/Users/mhall44/Mark/ibnw/OSPREY/dataFiles/rama500-prepro.data" } );
            
            //Create StructGrafter for structures struct1, struct2
            StructGrafter sg = new StructGrafter("../input_structs/4LSS.renum.pdb",
                    "../input_structs/4JM2.renum.pdb");
            
            //align struct2 to struct1
            sg.SSAlignment(new String[][] {new String[]{"791","797"},new String[]{"805","812"}}, 
                    new String[][] {new String[]{"91","97"},new String[]{"105","112"}},
                    "4LSS.4JM2.ssalign.pdb");
            //OR anchorAlignment
            //sg.anchorAlignment("797","802","97","102","4LSS.4JM2.anchoralign.pdb");
            
            //graft structure together
            //start with whole antibody part of 4LSS
            sg.initCompositeStructure("701","1616");
            
            //graft in loop
            sg.graftLoop("794", "804", "94", "104");//was 97-105...
            //add ligand
            sg.transferLigand("689","1127");
            
            //let's see what the structure is like now
            sg.struct3.connectivity12Valid = true;//let's play pretend first
            new KSParser().saveMolecule(sg.struct3, "4LSS.4JM2.open.pdb", 0);
            //These settings give a good 4LSS/4JM2 graft geometry- and BB clash-wise.  
            
            //close loops
            sg.matchLoopAnchors("794","800","800D","804","4LSS.4JM2.grafted.pdb");
            
            
            
            //make graft with CD4bs ligand ("graftc")
            StructGrafter sg2 = new StructGrafter("4LSS.4JM2.grafted.pdb","../input_structs/4LSS.renum.pdb");
            sg2.initCompositeStructure("701","1616");
            sg2.transferLigand("44","492");
            sg2.struct3.connectivity12Valid = true;
            new KSParser().saveMolecule(sg2.struct3, "4LSS.4JM2.graftc.pdb", 0);

            System.exit(0);*/
            
            /*EnvironmentVars.autoFix = false
            StructGrafter sg = new StructGrafter("4LSS.renum.pdb","4JM2.renum.pdb","4LSS.4JM2.graft3.pdb",
                    new String[][] {new String[]{"797","802"}},
                    new String[][] {new String[]{"97","102"}}, null);
            StructGrafter sg = new StructGrafter("4LSS.renum.pdb","4JAN.renum.pdb","4LSS.4JAN.graft3.pdb",
                    new String[][] {new String[]{"797","802"}},
                    new String[][] {new String[]{"797","802"}}, null);
            StructGrafter sg = new StructGrafter("4LSS.renum.pdb","3U2S.renum.pdb","4LSS.3U2S.graft3.pdb",
                    new String[][] {new String[]{"797","802"}},
                    new String[][] {new String[]{"997","1002"}}, null);*/
            
            //sg.graft();
            /*sg.SSAlignment(new String[][] {new String[]{"791","797"},new String[]{"805","812"}}, 
                    new String[][] {new String[]{"991","997"},new String[]{"1005","1012"}});//for 3U2S
            sg.SSAlignment(new String[][] {new String[]{"791","797"},new String[]{"805","812"}}, 
                    new String[][] {new String[]{"791","797"},new String[]{"805","812"}});//for 4JAN
            sg.SSAlignment(new String[][] {new String[]{"791","797"},new String[]{"805","812"}}, 
                    new String[][] {new String[]{"91","97"},new String[]{"105","112"}});//for 4JM2
                    */
            
            
            
            //System.exit(0);
            
            /*Object mat1 = KSParser.readObject("wt_rigidminM_COM.dat");
            Object mat2 = KSParser.readObject("wt_2str_rigidminM_COM.dat");
            
            System.exit(0);*/
            
            /*new VEGASTest().main();
            System.exit(0);*/
                    
            
            //DEBUG!!!
            /*SeqTreeNode seqNode = (SeqTreeNode)KSParser.readObject("BESTSEQNODE.dat");
            
            System.exit(0);*/
            
        
            //DEBUG!!!!
            /*NWPFSampleSet nss = (NWPFSampleSet)KSParser.readObject("NWPFSampleSetPROBLEM.dat");
            
            nss.updateFitVal(0);
            
            System.exit(0);*/
            
            
            
            /*MultivariatePoly mp = new MultivariatePoly(1,3);
            mp.terms.set(1.,0);
            mp.terms.set(2.);
            mp.times(mp).times(mp).print();
            mp.shift(DoubleFactory1D.dense.make(1,3)).print();
            mp.times(mp).times(mp).shift(DoubleFactory1D.dense.make(1,3)).print();
            
            
            MultivariatePoly mp2 = new MultivariatePoly(6,1);
            mp2.terms.set(2.5,3);
            mp2.terms.set(-0.5,2);
            
            System.out.println();
            mp2.print();
            DoubleMatrix1D c = DoubleFactory1D.dense.make(new double[] {0.1,0.2,0.3,0.4,0.5,0.6});
            mp2.shift(c).print();
            
            MultivariatePoly mp3 = new MultivariatePoly(6,1);
            mp3.terms.set(-2.,0);
            mp2.times(mp3).print();
            
            mp2.add(mp);
            mp2.print();
            
            
            MultivariatePoly mp2 = new MultivariatePoly(5,2);
            mp2.terms.set(1.,4,3);
            mp2.terms.set(1.,3,2);
            mp2.terms.set(1.,1);
            mp2.setConstant(7.);
            DoubleMatrix1D c = DoubleFactory1D.dense.make(new double[] {6,-2,0.5,3,0});
            System.out.println(c);
            mp2.print();
            mp2.shift(c).print();
            
            System.exit(0);*/
            
            
            
            
		KSParser parser = new KSParser();
		parser.checkMPI(args);
		System.out.println("KStar Finished.");

	}  // End main()
	
} // end class
