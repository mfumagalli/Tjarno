
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

You can create variables to store the path to these software.
For instance, these are mine:
```
MS=~/Software/msms/bin/msms
SAMTOOLS=samtools
ANGSD=~/Software/angsd
NGSTOOLS=~/Software/ngsTools
```
However you may no need to do that as when you load them they should be in your /usr/bin so you can call them directly without specifying the full path.

# Simulations with selection

We are using again the software [msms](http://www.mabs.at/ewing/msms/download.shtml) to perform coalescent simulations with selection.

Likewise, we are using again a model for the history of Europeans, East Asians and Native Americans following the studies [here](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695) and [here](http://www.ncbi.nlm.nih.gov/pubmed/26198033).

We are now simulating some genomic data (400kbp) for 3 populations, Europeans (EUR), East Asians (EAS) and Native Americans (NAM), with 10 diploid individuals each.

We are simulating that positive selection is acting on a derived allele of a variant roughly in the middle of this region.
This allele was not present in EUR, while at low frequency in EAS, and selection started after the split between EAS and EUR.

We can generate such data using msms with the following command:
```
$MS -N 7310 -ms 60 1 -t 300 -r 300 400000 -I 4 0 20 20 20 -n 1 1.68 -n 2 3.73 -n 3 7.29 -n 4 0.25 -eg 0 2 116 -eg 0 3 160 -ma x 0.88 0.56 0.00 0.88 x 2.79 0.00 0.56 2.79 x 0.00 0.00 0.00 0.00 x -ej 0.027 4 3 -ej 0.029 3 2 -en 0.029 2 0.29 -en 0.30 1 1 -SI 0.0137 4 0 0 0.01 0.01 -Sp 0.5 -Sc 0 1 0 0 0 -Sc 0 2 0 0 0 -Sc 0 3 0 0 0 -Sc 0 4 584 292 0-seed 1234 > Data/pops.ms
```
It's not really important how this line has been generated.
Look at the output file.
```
less -S Data/pops.ms
```




msms 3 pops, 500 kbp

# Site Frequency Spectra

single SFS

2D SFS

3D-SFS

# FST and PBS

FST and PBS

# Nucleotide diversity

TD 






