#Main directory: KStar.jar, Print.jar, read_prot.py, FastPD.exe, all of the protein files
#Protein files: pdb file, KStar.cfg, System.cfg DEE_pert.cfg

Syst=".pdbSystem.cfg"
DEE=".pdbDEE.cfg"

#for dir in 1AB2  1FSV  1NHZ  2EKO  2O2T; do 
for dir in 1AB2; do 
	cd $dir
#	java -jar ../KStar.jar -c KStar.cfg fixStruct $dir".pdb" $dir"_fixed.pdb"
        python ../checkPDB.py $dir".pdb" 3
	java -jar ../KStar.jar -c KStar.cfg doDEE ${dir}${Syst} ${dir}${DEE}
	for dat in `ls *minM*.dat`; do 
		java -jar ../Print.jar $dat
		python ../read_prot.py "min.dat"
		wine ../FastPD matrices.bin results.bin
	cd ../
	done
done
