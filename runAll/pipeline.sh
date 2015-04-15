#Main directory: KStar.jar, Print.jar, read_prot.py, FastPD.exe, all of the protein files
#Protein files: pdb file, KStar.cfg, System.cfg DEE_pert.cfg

Syst=".pdbSystem.cfg"
DEE=".pdbDEE.cfg"

#for dir in 1AB2  1FSV  1NHZ  2EKO  2O2T; do 
for dir in 1AB2; do 
	cd $dir
#	java -jar ../KStar.jar -c KStar.cfg fixStruct $dir".pdb" $dir"_fixed.pdb"
	echo "checkPDB start"
        python ../checkPDB.py $dir"_fixed.pdb" 3
	echo "KStar.jar start"
	java -jar ../KStar.jar -c KStar.cfg doDEE ${dir}${Syst} ${dir}${DEE}
	echo "Print.jar start"
	java -jar ../Print.jar min.dat
	echo "read_prot.py start"
	python ../read_prot.py min.dat_2d.txt
	echo "FastPD start"
	wine ../FastPD matrices.bin results.bin
	cd ../
done
