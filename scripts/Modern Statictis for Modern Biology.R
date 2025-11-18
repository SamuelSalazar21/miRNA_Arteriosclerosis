# tutorial de principales funciones estadisticas

HIV <- dpois(x=0:12, lambda = 5)
barplot(HIV, names.arg = 0:12, col="purple")

genotype=c("AA", "AO","BB", "A0", "OO", "AA", "BO", "BO","AO")
table(genotype)

genotypeF <- factor(genotype)
levels(genotypeF)
table(genotypeF)

help(factor)


rbinom(12, prob=2/3, size = 1)
rbinom(10, prob=0.3, size=15)
probability=dbinom(0:10^6, prob = 0.3, size=15)


round(probability,2)

barplot(probability)
1+1

