#Main directory: KStar.jar, Print.jar, read_prot.py, FastPD.exe, all of the protein files
#Protein files: pdb file, KStar.cfg, System.cfg DEE_pert.cfg

#Need to test!

for dir in 1AB2  1FSV  1NHZ  2EKO  2O2T; do 
#	echo $dir
	cd $dir
	java -jar KStar.jar -c KStar.cfg doDEE System.cfg DEE_pert.cfg
	for dat in `ls *minM*.dat`; do 
#	for dat in `ls *`; do
#		echo $dat"_2d"
		java -jar ../Print.jar $dat
		python ../read_prot.py $dat"_2d.txt"
		wine ../FastPD.exe matrices.bin results.bin
	cd ../
	done
done
