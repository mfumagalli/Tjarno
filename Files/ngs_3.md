
## Allele frequencies


**EXERCISE!!!**

Recalling our previous attempt with PCA, can we do even better?
Can you use the information on the estimated allele frequencies to ameliorate our inferences of population structure?
For instance, we can (i) use the estimates of allele frequencies to filter out sites which are likely to be monomorphic (hint: look at -minMaf option in ANGSD) and (ii) use the estimates of allele frequencies as prior information on our genotype probabilities (e.g. assuming HWE, hint: look at -doPost options).


solution on a different file

Use an informative prior and filter for minor allele frq
```
	$ANGSD/angsd -glf Data/pops.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doPost 1 -doGeno 32 -minMaf 0.02 -out Results/pops.inf.bin

	gunzip Results/pops.inf.bin.geno.gz

	NSITES=`zcat Results/pops.inf.bin.mafs.gz | tail -n +2 | wc -l`
	echo $NSITES
        
	$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/pops.inf.bin.geno -outfile Results/pops.inf.covar -nind $NIND -nsites $NSITES -call 0 -norm 0

        Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/pops.inf.covar -a Results/pops.clst -c 1-2 -o Results/pops.inf.pca.pdf

	evince Results/pops.inf.pca.pdf
```



