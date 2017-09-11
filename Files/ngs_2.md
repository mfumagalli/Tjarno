
## Population structure

Population structure
```
	$ANGSD/angsd -glf Data/pops.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doPost 2 -doGeno 32 -out Results/pops.flat.bin
        gunzip Results/pops.flat.bin.geno.gz
```

How many sites?
```
        NSITES=`zcat Results/pops.flat.bin.mafs.gz | tail -n +2 | wc -l`
        echo $NSITES
```

Standard PCA
```
        $NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/pops.flat.bin.geno -outfile Results/pops.flat.covar -nind $NIND -nsites $NSITES -call 1 -norm 1

        Rscript -e 'write.table(cbind(rep(seq(1,10),2),rep(seq(1,10),2),c(rep("EAS",10),rep("NAM",10))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="Results/pops.clst", quote=F)'

        Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/pops.flat.covar -a Results/pops.clst -c 1-2 -o Results/pops.flat.pca.pdf

```

Look at the plot:
```
	evince Results/pops.flat.pca.pdf
```

**QUESTION**
What's wrong here?

What can you do to improve your PCA?

```
	$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/pops.flat.bin.geno -outfile Results/pops.probs.covar -nind $NIND -nsites $NSITES -call 0 -norm 0

	Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/pops.probs.covar -a Results/pops.clst -c 1-2 -o Results/pops.probs.pca.pdf

	evince Results/pops.probs.pca.pdf
```

**QUESTION**

Can we do even better?

(filter based on geno probs and allele freq?)

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



## Additional

Genetic distances and plot a tree or do a MDS




