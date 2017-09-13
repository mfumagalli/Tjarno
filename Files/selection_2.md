
## Using ABC to quantify selection

Let's use a very simple approach to understand the possibiliy of using ABC to quantify selection parameters.
Assume that our goal is to estimate when selection started and how strong it was.
We want to use a simulation-based approach using the previously calculated statistics as summary statistics.
We will just use a very simple rejection-based method, as our goal is simply to illustrate the point of using ABC beyond demography inferences.

# Selection coefficient

For simplicity, let's first assume we know the selection time and that we want to estimate the selection coefficient in NAM.
We can use the FST and Tajima's D calculated on the 20kbp windows as summary statistics.

The first thing we need to decide is our prior distribution.
Imagine that we choose a uniform prior bounded at 0 (no selection) to 750, which, given our parametrisation on Ne, is roughly equal to a selection coefficient at the allele of 5%.

```
# define prior
MIN=0
MAX=750

# sample from the prior and calculate summary statistics
for i in `seq 1 100`;
do

  selcoeff=$RANDOM
  let "selcoeff %= $MAX"

  # sampled value
  echo $selcoeff
  selcoeffHomo=$(($selcoeff * 2)) # additive model

  # simulate data given this value for our parameter we want to estimate (you can see that -Sc is the flag to set the selection coefficient in the msms code)
  # note that we simulate 20kbp
  $MS -N 7310 -ms 40 1 -t 14 -r 14 20000 -I 4 0 0 20 20 -n 1 1.68 -n 2 3.73 -n 3 7.29 -n 4 0.25 -eg 0 2 116 -eg 0 3 160 -ma x 0.88 0.56 0.00 0.88 x 2.79 0.00 0.56 2.79 x 0.00 0.00 0.00 0.00 x -ej 0.027 4 3 -ej 0.029 3 2 -en 0.029 2 0.29 -en 0.30 1 1 -Sp 0.5 -SI 0.02 4 0 0 0.05 0.05 -Sc 0 3 0 0 0 -Sc 0 4 $selcoeffHomo $selcoeff 0 -seed 1234 > Data/tmp.ms
  
  # let's do something a bit different and instead of calculating summary statistics from known genotypes, let's do it from simulated NGS data, as this may lead to simulations closer to our real data
  DEPTH=4 # pick the same values you chose previously
  ERR=0.0075
  $ANGSD/misc/msToGlf -in Data/tmp.ms -out Data/tmp -regLen 20000 -singleOut 1 -depth $DEPTH -err $ERR -pileup 0 -Nsites 1 2> /dev/null
  $ANGSD/misc/splitgl Data/tmp.glf.gz 20 1 10 > Data/tmp1.glf.gz 2> /dev/null
  $ANGSD/misc/splitgl Data/tmp.glf.gz 20 11 20 > Data/tmp2.glf.gz 2> /dev/null



```




