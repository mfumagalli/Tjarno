
Open R.

Are the summary statistics correlated with the parameter?
```
values <- read.table("Results/simul.txt", head=T, stringsAsFactors=F)

cor(values)

pdf(file="Results/cor.pdf")
plot(values$SelCoeff, values$FST) 
plot(values$SelCoeff, values$TD_NAM)
dev.off()

```

Suppose we pick Tajima's D in NAM as the most informative summary statistics (it could be FST).
We can easily calculate the distance between each simulated value and our observation.
Still in R...
```
observed <- (-1.44)
dist <- abs(values$TD_NAM-observed)
# look at the distribution
hist(dist)
# select the top percentile (1st, 5th, 10th...)
topPerc <- 0.10*length(dist)
retained <- sort(dist, decreasing=F, index.return=T)$ix[1:topPerc]
# look at the distribution of simulated parameters corresponding to the retained summary statistics
pdf(file="Results/selcoeff.pdf")
hist(values$SelCoeff[retained])
# this is your posterior probability distribution for the selection coefficient
# you can even calculate the 90th interval
library(coda)
x <- values$SelCoeff[retained]
hpd <- HPDinterval(as.mcmc(x))
abline(v=hpd, lty=1)
dev.off()
# or plot a density
plot(density(x))
```




