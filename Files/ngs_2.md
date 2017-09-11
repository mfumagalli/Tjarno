
## Population structure

Now we know how to calculate genotype probabilities.
We want to use this information to assess the structure of our populations, for instance by doing a principal component analysis (PCA).
Again we can do that using ANGSD/ngsTools, either by calling or not calling genotypes.

The first thing we need to do is to calculate genotype posterior probabilities in binary mode using `-doGeno 32`:
```
$ANGSD/angsd -glf Data/pops.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doPost 2 -doGeno 32 -out Results/pops.flat.bin
gunzip Results/pops.flat.bin.geno.gz
```

**QUESTION** 
How many sites do we have?
```
NSITES=`zcat Results/pops.flat.bin.mafs.gz | tail -n +2 | wc -l`
echo $NSITES
```

The standard approach is to assign individual genotypes and calculate the covariance matrix as in Patterson et al. 2006.
This can be achieved by using ngsCovar package in ngsTools:
if you type:
```
$NGSTOOLS/ngsPopGen/ngsCovar
```
you will see a help message.

**QUESTION**
What information/files are required by reading the help message?
What's the command line to calculate suc covariance matrix?

-------------------------------------------------------------------------

Here is a possible solution.
You need to specify the number of sites and samples:
```
$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/pops.flat.bin.geno -outfile Results/pops.flat.covar -nind $NIND -nsites $NSITES -call 1 -norm 1
```

For plotting purposes, we need to create a plink cluster file to assign each sample to a population.
This is achieved using:
```
Rscript -e 'write.table(cbind(rep(seq(1,10),2),rep(seq(1,10),2),c(rep("EAS",10),rep("NAM",10))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="Results/pops.clst", quote=F)'
```

Finally, we are able to perform an eigenvector decomposition of the covariance matrix and plot the first two components:
```
Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/pops.flat.covar -a Results/pops.clst -c 1-2 -o Results/pops.flat.pca.pdf
```

**QUESTION**
Look at the plot:
```
evince Results/pops.flat.pca.pdf
```
What's wrong here?

--------------------------------------------------

**EXERCISE**
What can you do to improve your PCA?
Perhaps we can estimate the covariance matrix by not calling genotypes but rather using their probabilities.
Look at the help message at `$NGSTOOLS/ngsPopGen/ngsCovar` to see how you can achieve this.
Plot the new PCA and compare it with the previous one.
Try to plot other components (are they informative?).
A possible solution is available [here](https://github.com/mfumagalli/Tjarno/edit/master/Files/ngs_2_solution.md).




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



## Very optional: genetic distances

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





