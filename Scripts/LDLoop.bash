for i in {1..1000};
	do Rscript Scripts/LDSimCalcFunc.R 'Output/myseqdata'$i 40;
done
