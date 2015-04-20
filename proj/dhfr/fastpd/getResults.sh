for file in `ls *2d.txt`; do
	head -n 2 $file | tail -n 1  > "labels.txt"
        Rscript getLabels.R
	python ../../finalResults.py
done

