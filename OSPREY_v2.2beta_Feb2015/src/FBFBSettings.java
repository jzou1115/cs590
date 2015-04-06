
import java.util.ArrayList;
import java.util.BitSet;
import java.util.StringTokenizer;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author mhall44
 */
public class FBFBSettings {
    //This is to keep track of settings for FullBBFlexBlocks
    
    boolean useFBFB;
    
    double maxTrans;//how far do we allow translation in each dimension? (in angstroms)
    double maxRot;//max rotation for each axis (in degrees)
    double maxDihChange;//maximum change in each anchor dihedral (degrees)
    
    int numFBFB;//how many FullBBFlexBlocks there are
    
    //start and end residues for all the blocks.  These are PDB residue numbers
    String startRes[];
    String endRes[];
    
    boolean anchorPref;//anchor dihedral preferences (true=starting end) for all the blocks
    
    
    
    public FBFBSettings(){
        //by default, no FBFB
        useFBFB = false;
    }

    
    public FBFBSettings(ParamSet params){
        //initialize from input parameter set
        String FBFBRes = (String)params.getValue("FULLBBFLEX", "None");
        
        if(FBFBRes.equalsIgnoreCase("None"))//no FBFB
            useFBFB = false;
        else {
            //FBFBRes should be of the form "s1 e1 s2 e2..." where the starting and ending residues of
            //the FullBBFlexBlocks are (s1 to e1), (s2 to e2), ...
            
            useFBFB = true;
            
            StringTokenizer st = new StringTokenizer(FBFBRes);
            int numTokens = st.countTokens();
            if(numTokens%2==1)
                throw new RuntimeException("ERROR: Bad formatting of FULLBBFLEX");
            
            numFBFB = numTokens/2;
            startRes = new String[numFBFB];
            endRes = new String[numFBFB];
            
            int count = 0;
            while(st.hasMoreTokens()){
                startRes[count] = st.nextToken();
                endRes[count] = st.nextToken();
                count++;
            }
            
            maxTrans = (new Double((String)params.getValue("FBFBMAXTRANS","0.5"))).doubleValue();
            maxRot = (new Double((String)params.getValue("FBFBMAXROT","5"))).doubleValue();
            maxDihChange = (new Double((String)params.getValue("FBFBMAXDIHCHANGE","10"))).doubleValue();
            
            
            anchorPref = (new Boolean((String)params.getValue("FBFBANCHORPREF", "false"))).booleanValue();
            
            boolean useCCD = (new Boolean((String)params.getValue("USECCD", "false"))).booleanValue();
            if(!useCCD)
                throw new RuntimeException("ERROR: FBFB not supported without CCD.");
        }  
    }
    
    
    
    public FullBBFlexBlock[] makeFBFB(Molecule m){
        //create any of the FullBBFlexBlocks that are applicable to residues in m
        //we'll assume that if both the start and end residues are present, then the block is applicable
        //(unapplicability expected to occur for unbound subunits during K*, etc.)
        
        ArrayList<FullBBFlexBlock> fList = new ArrayList<>();
        
        for(int f=0; f<numFBFB; f++){
            int blockStartRes = m.mapPDBresNumToMolResNum(startRes[f]);
            int blockEndRes = m.mapPDBresNumToMolResNum(endRes[f]);
            
            if(blockStartRes>-1 && blockEndRes>-1)
                fList.add( new FullBBFlexBlock(m,blockStartRes,blockEndRes,anchorPref,maxTrans) );
            //we'll restrict closure block CA motions using maxTrans (which of course also regulates peptide plane translations)
        }
        
        return fList.toArray(new FullBBFlexBlock[fList.size()]);
    }
}
