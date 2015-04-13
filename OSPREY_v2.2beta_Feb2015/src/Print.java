import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class Print{	

	public static void main(String[] args){

		try{
			//String file= "/home/jennifer/Documents/2015Spring/CS590/finalProject/OSPREY_v2.2beta_Feb2015/example/dhfr-ks-example/dhfr-dee-examplemaxM_0.dat";
			String file= args[0];
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file));

           // PairwiseEnergyMatrix arpMatrix = new PairwiseEnergyMatrix();
            //arpMatrix.eMatrix = (double [][][][][][])in.readObject();
            double[][][][][][] emat = (double [][][][][][])in.readObject();
            //double[][][][][][] emat= arpMatrix.eMatrix;
            
			in.close();
			
			if (emat==null){
				System.out.println("Error: eMatrix null; file not written");
				return;
			}
			
			
			BufferedWriter fw2= new BufferedWriter(new FileWriter(file+"_deepToString.txt"));
			fw2.write(Arrays.deepToString(emat));
			fw2.close();
		
		
			HashMap<String, Integer> labelsMap= new HashMap<String, Integer>(); 
			List<String> labels = new ArrayList<String>();
			int i=0;
			for (int p1=0; p1<emat.length-1; p1++){
				if (emat[p1]!=null){
					for (int a1=0; a1<emat[p1].length; a1++){
						if (emat[p1][a1]!=null){
							for (int r1=0; r1<emat[p1][a1].length; r1++){
								if (emat[p1][a1][r1]!=null){
									String line="";
									line= p1+"_"+a1+"_"+r1;
									if(!labelsMap.keySet().contains(line)){
										labelsMap.put(line, i);
										labels.add(line);
										i++;
									}
									else{
										System.out.println("Error: Duplicate Key "+ line);
									}
									/**
									for(int p2=0; p2<emat[p1][a1][r1].length;p2++){
										if(emat[p1][a1][r1][p2]!=null){
											for(int a2=0; a2<emat[p1][a1][r1][p2].length; a2++){
												if(emat[p1][a1][r1][p2][a2]!=null){
													for(int r2=0; r2<emat[p1][a1][r1][p2][a2].length;r2++){
														String line2="";
														line2= p2+"_"+a2+"_"+r2;
														if(!labelsMap.keySet().contains(line2)){
															labelsMap.put(line2, i);
															labels.add(line2);
															i++;
														}
														else{
															System.out.println("Duplicate Key (not error)"+ line);
														}
													}
												}
											}
										}
									}
									**/
								}
								
							}
						}
					}
				}				
			}
			//System.out.println("SIZE OF EMATRIX IS: "+labels.size());
			
			double[][] out= new double[labels.size()][labels.size()];
		

			for(String l:labels){
				for(String k:labels){
					String[] tokens1= l.split("_");
					String[] tokens2= k.split("_");
					int p1= Integer.parseInt(tokens1[0]);
					int a1= Integer.parseInt(tokens1[1]);
					int r1= Integer.parseInt(tokens1[2]);
					int p2= Integer.parseInt(tokens2[0]);
					int a2= Integer.parseInt(tokens2[1]);
					int r2= Integer.parseInt(tokens2[2]);
					if(!labelsMap.containsKey(l) | !labelsMap.containsKey(k)){
						System.out.println("Missing key!\t"+tokens1+"\t"+tokens2);
					}
					else if(p1==p2){
						try{
						if(r1==r2 && a1==a2){
							out[labelsMap.get(l)][labelsMap.get(k)]= emat[p1][a1][r1][p1][0][0]+emat[p1][a1][r1][p1][0][1];	
						}
						else{
						out[labelsMap.get(l)][labelsMap.get(k)]= 0;
						}
						}catch (Exception e){
							e.printStackTrace();
						}
					}
					else{
						//null pointer exception at third conditional, even though second passed
						
						if((emat!=null) && (emat.length>p1)&&(emat[p1]!=null) && (emat[p1].length>a1) && (emat[p1][a1]!=null) &&(emat[p1][a1].length>r1) && (emat[p1][a1][r1]!=null) && (emat[p1][a1][r1].length>p2) && (emat[p1][a1][r1][p2]!=null) && (emat[p1][a1][r1][p2].length>a2) && (emat[p1][a1][r1][p2][a2]!=null) && (emat[p1][a1][r1][p2][a2].length>r2)){
							//System.out.println("e="+emat[p1][a1][r1][p2][a2][r2]);
							out[labelsMap.get(k)][labelsMap.get(l)] = emat[p1][a1][r1][p2][a2][r2];
							out[labelsMap.get(l)][labelsMap.get(k)] = emat[p1][a1][r1][p2][a2][r2];
						}
					}
				}
			}
			
			BufferedWriter fw= new BufferedWriter(new FileWriter(file+"_2d.txt"));

			List<Integer> rotCounts= new ArrayList<Integer>();
			int r=0;
			int p=0;
			String lab="";
			for(String l:labels){
				int pos= Integer.parseInt(l.split("_")[0]);
				if(pos ==p){
					r++;
				}
				else{
					rotCounts.add(r);
					p++;
					r=1;
				}
				lab=lab+l+"\t";
			}
			rotCounts.add(r);
			String rc="";
			for(int a=0; a<rotCounts.size();a++){
				rc=rc+rotCounts.get(a)+"\t";
			}
			fw.write(rc.trim()+"\n");
			fw.write(lab.trim()+"\n");
			
			int row=0;
			for(String l: labels){
				String line= "";
				for(int k=0; k<labels.size();k++){
					line=line+out[row][k]+"\t";
					//System.out.println(line);
				}
				line=line.trim()+"\n";
				fw.write(line);
				row++;
			}
			fw.close();
			
			return;
		}
		catch (Exception e){
			e.printStackTrace();
		}		
	}
}
