
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




