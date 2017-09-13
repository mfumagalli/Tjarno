
## Inference of selection from NGS data

# Preparation

Create some folders where to put stuff.
```
mkdir Selection
cd Selection
mkdir Data
mkdir Results
```

Remember to load all required software.
```
module load msms/v3.2rc-b163
module load samtools/v1.3.1
module load angsd/v0.918
module load ngsTools/vX.X
```
Make sure you use ANGSD version 0.918.

You can create variables to store the path to these software.
For instance, these are mine:
```
MS=~/Software/msms/bin/msms
SAMTOOLS=~/Software/samtools-1.5/samtools
ANGSD=~/Software/angsd
NGSTOOLS=~/Software/ngsTools
```
However you may no need to do that as when you load them they should be in your /usr/bin so you can call them directly without specifying the full path.

# Simulations with selection

We are using again the software [msms](http://www.mabs.at/ewing/msms/download.shtml) to perform coalescent simulations with selection.

Likewise, we are using again a model for the history of East Asians and Native Americans following the studies [here](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695) and [here](http://www.ncbi.nlm.nih.gov/pubmed/26198033).

We are now simulating some genomic data (500kbp) for 2 populations, East Asians (EAS) and Native Americans (NAM), with 10 diploid individuals each.

We are simulating that positive selection is acting on a derived allele of a variant roughly in the middle of this region.
This allele was not present in EUR, while at low frequency in EAS, and selection started after the split between EAS and EUR.

We can generate such data using msms with the following command:
```
$MS -N 7310 -ms 40 1 -t 350 -r 350 500000 -I 4 0 0 20 20 -n 1 1.68 -n 2 3.73 -n 3 7.29 -n 4 0.25 -eg 0 2 116 -eg 0 3 160 -ma x 0.88 0.56 0.00 0.88 x 2.79 0.00 0.56 2.79 x 0.00 0.00 0.00 0.00 x -ej 0.027 4 3 -ej 0.029 3 2 -en 0.029 2 0.29 -en 0.30 1 1 -Sp 0.5 -SI 0.02 4 0 0 0.05 0.05 -Sc 0 3 0 0 0 -Sc 0 4 1000 500 0 -seed 1234 > Data/pops.ms
```
It's not really important how this line has been generated.
Look at the output file.
```
less -S Data/pops.ms
```

If you recall what did the other day, we have to generate and index a "fake" reference ancestral sequence:
```
Rscript -e 'cat(">ref\n",paste(rep("A",5e5),sep="", collapse=""),"\n",sep="")' > Data/ref.fa 
$SAMTOOLS faidx Data/ref.fa
```

Again, we now  simulate genomic data by setting the average depth per site per sample (pick sometyhing between 2 and 8) and the sequencing error rate (between 0.1% and 1%).
We use a utility in ANGSD to simulate sequencing reads based on these parameters and split the genotype likelihoods into the 3 populations:
```
DEPTH=4 # pick a value between 2 and 8
ERR=0.0075

$ANGSD/misc/msToGlf -in Data/pops.ms -out Data/pops -regLen 500000 -singleOut 1 -depth $DEPTH -err $ERR -pileup 0 -Nsites 1 -seed 1234

$ANGSD/misc/splitgl Data/pops.glf.gz 20 1 10 > Data/pop1.glf.gz 
$ANGSD/misc/splitgl Data/pops.glf.gz 20 11 20 > Data/pop2.glf.gz 

```

# Site Frequency Spectra

As we discussed the other day, a good prior for estimating the per-site allele frequencies is the Site Frequency Spectrum (SFS).
This is basis for deriving any summary statistics based on allele frequencies.

Let's estimate the SFS byt first calculating the sample allele frequency likelihoods in the .saf files.
We use the utility `realSFS` to estimate the SFS.
```
NIND=10
for i in `seq 1 2`; do $ANGSD/angsd -glf Data/pop$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doSaf 1 -out Results/pop$i; $ANGSD/misc/realSFS Results/pop$i.saf.idx > Results/pop$i.sfs; done
```

**QUESTION**
How many values do you expect?

If you open one file:
```
cat Results/pop1.sfs
```
you will see the estimated proportion of sites with a certain derived allele frequency.

We can parse this file to plot the SFS for each population:
```
Rscript $NGSTOOLS/Scripts/plotSFS.R Results/pop1.sfs-Results/pop2.sfs EAS-NAM 0 Results/pops.sfs.pdf
evince Results/pops.sfs.pdf
```

**QUESTION** Do they behave as expected (given what we know about human evolution)?

## multi-SFS

We now want to estimate a multi-dimensional SFS, for instance the joint SFS between 2 populations (2D). 
This can be used for making inferences on their divergence process (split time, migration rate and so on). 
We will use this information as prior probability for our FST estimation.

```
$ANGSD/misc/realSFS Results/pop1.saf.idx Results/pop2.saf.idx > Results/pops.2dsfs
```

**QUESTIONS**
How do the resulting files look like?
```
less -S Results/pops.2dsfs
```
The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.

You can plot it, but you need to define how many samples (individuals) you have per population.
```
Rscript $NGSTOOLS/Scripts/plot2DSFS.R Results/pops.2dsfs EAS-NAM 10-10
evince Results/pops.2dsfs.pdf
```

## FST

A nice summary statistics to quantify the 2D-SFS is FST, index of population genetic differentiation.
The local FST is highly informative in selection patterns.
```
$ANGSD/misc/realSFS fst index Results/pop1.saf.idx Results/pop2.saf.idx -sfs Results/pops.2dsfs -fstout Results/pops -whichFst 1

$ANGSD/misc/realSFS fst stats Results/pops.fst.idx
```
The weighted value is the estimate of FST for the whole region.
Take a note of this value.

We can do a sliding windows scan:
```
$ANGSD/misc/realSFS fst stats2 Results/pops.fst.idx -win 20000 -step 5000 > Results/pops.fst.wins
```
If you open the file you can see whether there are some windows have higher values.
```
less -S Results/pops.fst.wins
```
Look for the peak and record the highest value.

# Nucleotide diversity

To confirm the selection is acting in NAM samples, we can calculate some measure of genetic diversity, as Tajima's D.

First we calculate Tajima's D along windows.
```
NIND=10

for i in 1 2; do $ANGSD/angsd -glf Data/pop$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doSaf 1 -doThetas 1 -pest Results/pop$i.sfs -out Results/pop$i; $ANGSD/misc/thetaStat do_stat Results/pop$i.thetas.idx -win 20000 -step 5000; done
```

Have a look at Tajima's D values for both populations, for instance NAM, indexed as 2, at column 9.
```
cut Results/pop1.thetas.idx.pestPG -f 3,9 | less -S
```
Take a note of the lowest value.

**QUESTIONS** Do you observe any particular pattern?

Once we identify potential targets of selection, we may want to test for it and quantify its features.
We will use these summary statistics for this purpose.

Before moving on, let's compute Tajima's D for the whole region.
```
for i in 1 2; do $ANGSD/angsd -glf Data/pop$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doSaf 1 -doThetas 1 -pest Results/pop$i.sfs -out Results/pop$i; $ANGSD/misc/thetaStat do_stat Results/pop$i.thetas.idx; done
```
Look at Tajima's D values:
```
cut Results/pop?.thetas.idx.pestPG -f 9
```

**QUESTION** Is that expected? What's happening here?

