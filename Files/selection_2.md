
## Using ABC to quantify selection

Let's use a very simple approach to understand the possibiliy of using ABC to quantify selection parameters.
Assume that our goal is to estimate when selection started and how strong it was.
We want to use a simulation-based approach using the previously calculated statistics as summary statistics.
We will just use a very simple rejection-based method, as our goal is simply to illustrate the point of using ABC beyond demography inferences.
However you can use the same framework by doing more complex ABC inferences.

# Selection coefficient

For simplicity, let's first assume we know the selection time and that we want to estimate the selection coefficient in NAM.
We can use the FST and Tajima's D calculated on the 20kbp windows as summary statistics.

Let's initialise an empty file to store the results
```
echo SelCoeff FST TD_EAS TD_NAM > Results/simul.txt
```
where each row will have a value for our parameter and the corresponding summary statistics.

**1) Prior distribution**

The first thing we need to decide is our prior distribution.
Imagine that we choose a uniform prior bounded at 0 (no selection) to 300, which, given our parametrisation of selection coefficient equal to 2 times Ne with Ne=7310, is roughly equal to a selection coefficient at the allele of 2%.
However you can choose whatever range you want.
We also need to choose how many samplings to do,
Ideally, you need a lot of them but for time constrains, let's just do 100.
```
# 1) define prior, uniform
MIN=0
MAX=300
# 1b) how many simulations?
NSIM=100
```

**2) Simulate (sample from the prior) and calculate summary statistics**

```
for s in `seq 1 $NSIM`;
do

  # 2a) sampling from the prior
  # this is a trick to get bounded pseudo-random numbers in bash
  selcoeff=$RANDOM
  let "selcoeff %= $MAX"
  selcoeffHomo=$(($selcoeff * 2)) # additive model, this is the selection coefficient at the homozygous state

  # 2b) simulate data given this value for our parameter we want to estimate (you can see that -Sc is the flag to set the selection coefficient in the msms code)
  # note that we simulate 20kbp
  $MS -N 7310 -ms 40 1 -t 14 -r 14 20000 -I 4 0 0 20 20 -n 1 1.68 -n 2 3.73 -n 3 7.29 -n 4 0.25 -eg 0 2 116 -eg 0 3 160 -ma x 0.88 0.56 0.00 0.88 x 2.79 0.00 0.56 2.79 x 0.00 0.00 0.00 0.00 x -ej 0.027 4 3 -ej 0.029 3 2 -en 0.029 2 0.29 -en 0.30 1 1 -Sp 0.5 -SI 0.02 4 0 0 0.01 0.05 -Sc 0 3 0 0 0 -Sc 0 4 $selcoeffHomo $selcoeff 0 -seed 1234 > Data/tmp.ms
  
  # 2c) calculate simulated summary statistics
  
  # let's do something a bit different and instead of calculating summary statistics from known genotypes, let's do it from simulated NGS data, as this may lead to simulations closer to our real data (nor more noise...)
  DEPTH=4 # pick the same values you chose previously
  ERR=0.0075
  $ANGSD/misc/msToGlf -in Data/tmp.ms -out Data/tmp -regLen 20000 -singleOut 1 -depth $DEPTH -err $ERR -pileup 0 -Nsites 1 2> /dev/null
  $ANGSD/misc/splitgl Data/tmp.glf.gz 20 1 10 > Data/tmp1.glf.gz 2> /dev/null
  $ANGSD/misc/splitgl Data/tmp.glf.gz 20 11 20 > Data/tmp2.glf.gz 2> /dev/null

  # these are the commands already discussed to calculate FST and Tajima's D

  # sfs
  NIND=10
  for i in `seq 1 2`; do $ANGSD/angsd -glf Data/tmp$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doMajorMinor 4 -doMaf 1 -doSaf 1 -out Results/tmp$i 2> /dev/null; $ANGSD/misc/realSFS Results/tmp$i.saf.idx 2> /dev/null > Results/tmp$i.sfs; done

  $ANGSD/misc/realSFS Results/tmp1.saf.idx Results/tmp2.saf.idx 2> /dev/null > Results/tmp.2dsfs

  # fst
  $ANGSD/misc/realSFS fst index Results/tmp1.saf.idx Results/tmp2.saf.idx -sfs Results/tmp.2dsfs -fstout Results/tmp -whichFst 1 2> /dev/null
  FST=`$ANGSD/misc/realSFS fst stats Results/tmp.fst.idx 2> /dev/null | cut -f 2 -d " "`

  # Tajima's D
  TD=(NA NA)
  for i in 1 2; do $ANGSD/angsd -glf Data/tmp$i.glf.gz -ref Data/ref.fa -fai Data/ref.fa.fai -isSim 1 -nInd $NIND -doSaf 1 -doThetas 1 -pest Results/tmp$i.sfs -out Results/tmp$i 2> /dev/null; $ANGSD/misc/thetaStat do_stat Results/tmp$i.thetas.idx 2> /dev/null; TD[$((i-1))]=`cut Results/tmp$i.thetas.idx.pestPG -f 9 | tail -n 1`; done

  # now we have all we need, our sampled parameter and the corresponding summary statistics
  # let's save these values on a file
  echo $selcoeff $FST ${TD[*]} >> Results/simul.txt

  if [ $s -eq 1 ];
  then
    echo SIM COEFF FST TD
  fi
  
  echo $s $selcoeff $FST ${TD[1]}

done

```

If you are bored to wait, I provide you with a precompiled file [here](https://github.com/mfumagalli/Tjarno/edit/master/Files/simul.txt).
If by clicking on it doesn't work, you can clone the whole repository and move the file into your folder:
```
git clone https://github.com/mfumagalli/Tjarno
cp Tjarno/Files/simul.1K.txt Results/.
```
In this file I have 1,000 simulations.

Now you have everything you need to do your inferences.
In the simplest scenario, you can calculate the Eucledian distance between your observed summary and simulated statistics.
Then you retain the simulated parameters with the lowest distance.
Their distribution is your approximation of the posterior distribution of your parameter, the selection coefficient.

**EXERCISE** Using R or even a spreadsheet:
- assess the correlation between each summary statistic and the parameter
- select the most informative summary statistic
- use the rejection method to derive a (very much) approximated distribution of the selection coefficient (calculate distances between observed and simulated summary statistics, selected the smallest distances and plot the corresponding selection coefficient values)\\

A possible quick solution can be found [here](https://github.com/mfumagalli/Tjarno/edit/master/Files/selection_2_solution.md).

Comment on the results (e.g. is the prior appropriate? do we have power to estimate the parameter? does the posterior distribution follow our prior distribution? can we gain more power if we use multiple summary statistics?...)

# Selection coefficient and time

The joint estimation of selection strength and time is challenging, as these two parameters are tighly connected.
For instance, very high local values of FST may be compatible to both strong+recent selection and weak+old selection.
The key is to find a set of summary statistics that it's informative on the different behaviour of these selective scenatios.

**QUESTION** Which summary statistics would you use to distinguish between old and recent selection? Think of the effect of mutation and recombination on genetic diversity.


