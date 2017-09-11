
## Allele frequencies




## Sample allele frequencies likelihoods and probabilities

For each population, we are now calculating the sample allele frequency probabilities (what we call .saf files).
```
	NIND=10
	for i in `seq 1 2`;
        do
        	$ANGSD/angsd -glf Data/pop$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doSaf 1 -out Results/pop$i
	done
```

Let's inspect the files (for instance for population 1):
```
	$ANGSD/misc/realSFS print Results/pop1.saf.idx | less -S
```

**QUESTION**
What are these values?
Can you identify any SNP?

----------------------------------------------------------------------------

## (brief) SNP calling

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


------------------------------------------------------------------------------



Now we estimate some summary statistics.
For instance, we can calculate the expected number of polymorphic sites in our region, in both populations separately.
For doing that, let's calculate the probability of each site being variable (for instance, for the first 100 sites).
```
	NSITES=100
	NIND=10

	for i in `seq 1 2`;
	do
		zcat Results/pop$i.saf.gz > Results/pop$i.saf

		$NGSTOOLS/ngsPopGen/ngsStat -npop 1 -postfiles Results/pop$i.saf -outfile Results/pop$i.stats -nind $NIND -nsites $NSITES -iswin 0

	done
```

**QUESTION**

Are these values consistent with the sample allele frequency likelihoods?

What's annoying here?

If you calculate the expected number of polymorphic sites (use `-iswin 1 -block_size $NSITES`), do you obtain sensible value?
```
	i=1
	$NGSTOOLS/ngsPopGen/ngsStat -npop 1 -postfiles Results/pop$i.saf -outfile Results/pop$i.whole.stats -nind $NIND -nsites $NSITES -iswin 1 -block_size $NSITES
```

How can we solve this?


--------------------------------------------------------------------

We can estimate the site frequency spectrum (SFS) to be used as prior information.
Let's calculate the SFS for each population separately.

```
	for i in `seq 1 2`;
       	do
		$ANGSD/misc/realSFS Results/pop$i.saf.idx > Results/pop$i.sfs
	done
```

**QUESTION**

How many values do you expect?


Open the files:
```
	i=1
	cat Results/pop$i.sfs
```
What do these values represent?

We can plot them:
```
	Rscript $NGSTOOLS/Scripts/plotSFS.R Results/pop1.sfs-Results/pop2.sfs EAS-NAM 0 Results/pops.sfs.pdf
	evince Results/pops.sfs.pdf
```

**QUESTION**

Do they behave as expected (given what we know about human evolution)?


----------------------------------------------------------------------

Now we can use such SFS as prior information to improve our estimation of genetic diversity.
Let's calculate some diversity indexes.

```
	NIND=10
	for i in `seq 1 2`;
	do
		$ANGSD/angsd -glf Data/pop$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doSaf 1 -doThetas 1 -pest Results/pop$i.sfs -out Results/pop$i
		$ANGSD/angsd -bam bam.filelis -out out -doThetas 1 -doSaf 1 -pest out.sfs -anc chimpHg19.fa -GL
		$ANGSD/misc/thetaStat do_stat Results/pop$i.thetas.idx
	done
```

Open
```
	cat Results/pop?.thetas.idx.pestPG
```

**QUESTION**

Comment on theta and Tajima's D.

------------------------------------------------------------------------------------

**EXERCISE**

Do a sliding window for Tajima's D and plot it.

give tips and point to appropriate manual

then give link to a possible solution


-----------------------------------------------------------------------------------

## Additional

LD decay


