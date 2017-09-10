
## Preparation

Please make sure to follow these preparatory instructions:
```
mkdir Data
mkdir Results
```
```
MS=~/Software/msms/bin/msms
SAMTOOLS=~/Software/samtools-1.5/samtools
ANGSD=~/Software/angsd
NGSTOOLS=~/Software/ngsTools
```

We are using the software [msms](http://www.mabs.at/ewing/msms/download.shtml) to perform coalescent simulations under neutrality (and with selection, as we will see later).
Please follow the link to get the manual, if interested (but not required).

We also use a model previously estimated [here](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695) for the evolution of Africans, Europeans and East Asians.

We are adding the history of Native Americans to this model, roughly following estimates reported in [this](http://www.ncbi.nlm.nih.gov/pubmed/26198033) paper.
Thus, we will assume that Native Americans (their ancestors) splitted from East Asians 20kya and their effective population size is 2,000 from the split until present.

We are now simulating some genomic data (2Mbp) for 2 populations of 10 diploid individuals each.

We can generate such data with the following command:
```
	$MS -ms 40 1 -t 1500 -r 1500 2000000 -I 4 0 0 20 20 -n 1 1.68 -n 2 3.73 -n 3 7.29 -n 4 0.25 -eg 0 2 116 -eg 0 3 160 -ma x 0.88 0.56 0.00 0.88 x 2.79 0.00 0.56 2.79 x 0.00 0.00 0.00 0.00 x -ej 0.027 4 3 -ej 0.029 3 2 -en 0.029 2 0.29 -en 0.30 1 1 -seed 1234 > Data/pops.ms
```

Let's produce a reference sequence (with all As) and index it:
```
	Rscript -e 'cat(">ref\n",paste(rep("A",2e6),sep="", collapse=""),"\n",sep="")' > Data/ref.fa 
	$SAMTOOLS faidx Data/ref.fa
```

Now we need to simulate sequencing data.
You can decide the average depth per site per sample (pick sometyhing between 2 and 8) and the sequencing error rate.
```
	DEPTH=4
	ERR=0.0075
	$ANGSD/misc/msToGlf -in Data/pops.ms -out Data/pops -regLen 2000000 -singleOut 1 -depth $DEPTH -err $ERR -pileup 0 -Nsites 0 -seed 1234
```
and then we have to split the results into two populations:
```
	$ANGSD/misc/splitgl Data/pops.glf.gz 20 1 10 > Data/pop1.glf.gz 
	$ANGSD/misc/splitgl Data/pops.glf.gz 20 11 20 > Data/pop2.glf.gz 
```

The most intuitive way to look at these genotype likelihoods is to calculate posterior probabilities using a flat uniform prior.
We will focus only on the 3 possible genotypes after having inferred the two most likely alleles.

```
	NIND=20
        $ANGSD/angsd -glf Data/pops.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doPost 2 -doGeno 13 -out Results/pops.flat
```

Look for low probabilities or missing data.
```
        less -S Results/pops.flat.geno.gz
```

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



	










