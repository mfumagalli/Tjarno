
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

-------------------------------------------------------------

## Summary statistics

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





## 2D-SFS

We now want to estimate a multi-dimensional SFS, for instance the joint SFS between 2 populations (2D). 
This can be used for making inferences on their divergence process (split time, migration rate and so on). 
Here we are interested in estimating the 2D-SFS as prior information for our FST estimation.

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.

```
	$ANGSD/misc/realSFS Results/pop1.saf.idx Results/pop2.saf.idx > Results/pops.2dsfs
```

**QUESTIONS**

```
	less -S Results/pops.2dsfs
```

What do the values represent?

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.

You can plot it, but you need to define how many samples (individuals) you have per population.

```
	Rscript $NGSTOOLS/Scripts/plot2DSFS.R Results/pops.2dsfs EAS-NAM 10-10
	evince Results/pops.2dsfs.pdf
```

**QUESTION**

Comment.

What would you get if you had two populations very close to each other?
```
	$ANGSD/misc/realSFS Results/pop1.saf.idx Results/pop1.saf.idx > Results/fake.2dsfs
	Rscript $NGSTOOLS/Scripts/plot2DSFS.R Results/fake.2dsfs EAS-EAS 10-10
        evince Results/fake.2dsfs.pdf
```

----------------------------------------------------------------------`

A nice summary statistics to quantify the 2D-SFS is FST.

```
	$ANGSD/misc/realSFS fst index Results/pop1.saf.idx Results/pop2.saf.idx -sfs Results/pops.2dsfs -fstout Results/pops

```

```
	$ANGSD/misc/realSFS fst print Results/pops.fst.idx
```

```
	$ANGSD/misc/realSFS fst stats Results/pops.fst.idx
```


**EXERCISE**

Do a sliding window a plot.




