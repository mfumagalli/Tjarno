
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



## Additional: genetic distances

We can compute genetic distances as a basis for population clustering driectly from genotype probabilities, and not from assigned genotypes as we have seen how problematic these latters can be at low-depth.

First, we compute genotype posterior probabilities jointly for all samples using `-doGeno 8`:
```
$ANGSD/angsd -glf Data/pops.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doPost 2 -doGeno 8 -out Results/pops.flat.bin
```

Next we record how many sites we retrieve.
```
NSITES=`zcat Results/pops.flat.bin.mafs.gz | tail -n +2 | wc -l`
echo $NSITES
```

Then we create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("EAS","NAM"),each=10), rep(1:10, 2), sep="_"), sep="\n", file="Data/pops.label")'
cat Data/pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
```
$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/pops.flat.bin.geno.gz -probs -n_ind 20 -n_sites $NSITES -labels Data/pops.label -o Results/pops.dist

less -S Results/pops.dist
```

We can visualise the pairwise genetic distances in form of a tree.
For doing so you need to install [FastMe](http://www.atgc-montpellier.fr/fastme/).
```
FASTME=~/Software/fastme-2.1.5-linux64
$FASTME -D 1 -i Results/pops.dist -o Results/pops.tree -m b -n b
cat Results/pops.tree
```
Finally, we plot the tree.
```
Rscript $NGSTOOLS/Scripts/plotTree.R Results/pops.tree
evince Results/pops.tree.pdf
```





