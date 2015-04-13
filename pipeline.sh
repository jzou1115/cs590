for dir in protein1 protein2; do 
	cd $dir
	python read_prot.py
	wine FastPD.exe matrices.bin results.bin
	cd ../
done
