library(ggplot2)

#######
data <- read.table("cell_growth_rates_by_temp - concat4R.tsv", header=T)

########
#pdf("cells_per_day.pdf")
#ggplot(data, aes(x=cell_line, y=cells.day)) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#dev.off()

#pdf("cells_per_day.log10.pdf")
#ggplot(data, aes(x=cell_line, y=log10(cells.day))) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#dev.off()

#########

pdf("cell_doubling_time.pdf")
ggplot(data, aes(x=cell_line, y=doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("cell_doubling_time.ylim.pdf")
ggplot(data, aes(x=cell_line, y=doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,7)
dev.off()

pdf("cell_doubling_time.log10.pdf")
ggplot(data, aes(x=cell_line, y=log10(doubling_time_Td..t2.t1..ln.2..ln.q2.q1.))) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##########

pdf("cell_growth_per_day.pdf")
ggplot(data, aes(x=cell_line, y=growth_rate_per_day)) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("cell_growth_per_day.log10.pdf")
ggplot(data, aes(x=cell_line, y=log10(growth_rate_per_day))) + geom_boxplot() + geom_point(col=data$temp_..C., show.legend=TRUE, position="jitter") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##########

JW18DOX <- subset(data, cell_line=="JW18DOX")
JW18wMel <- subset(data, cell_line=="JW18wMel")
JW18wRi <- subset(data, cell_line=="JW18wRi")
JW18wWil <- subset(data, cell_line=="JW18wWil")

S2DOX <- subset(data, cell_line=="S2DOX")
S2wMel <- subset(data, cell_line=="S2wMel")
S2wRi <- subset(data, cell_line=="S2wRi")
S2wWil <- subset(data, cell_line=="S2wWil")

Dsim6B <- subset(data, cell_line=="Dsim6B")
Dsim6BwMel <- subset(data, cell_line=="Dsim6BwMel")
Dsim6BwWil <- subset(data, cell_line=="Dsim6BwWil")

###################
##### by temp #####
### 23C ###
JW18DOX23C <- subset(data, cell_line=="JW18DOX" & temp_..C.=="23")
JW18wMel23C <- subset(data, cell_line=="JW18wMel" & temp_..C.=="23")
JW18wRi23C <- subset(data, cell_line=="JW18wRi" & temp_..C.=="23")

S2DOX23C <- subset(data, cell_line=="S2DOX" & temp_..C.=="23")
S2wMel23C <- subset(data, cell_line=="S2wMel" & temp_..C.=="23")
S2wRi23C <- subset(data, cell_line=="S2wRi" & temp_..C.=="23")

Dsim6B23C <- subset(data, cell_line=="Dsim6B" & temp_..C.=="23")
Dsim6BwMel23C <- subset(data, cell_line=="Dsim6BwMel" & temp_..C.=="23")
Dsim6BwWil23C <- subset(data, cell_line=="Dsim6BwWil" & temp_..C.=="23")

### 26C ###
JW18DOX26C <- subset(data, cell_line=="JW18DOX" & temp_..C.=="26")
JW18wMel26C <- subset(data, cell_line=="JW18wMel" & temp_..C.=="26")
JW18wRi26C <- subset(data, cell_line=="JW18wRi" & temp_..C.=="26")

S2DOX26C <- subset(data, cell_line=="S2DOX" & temp_..C.=="26")
S2wMel26C <- subset(data, cell_line=="S2wMel" & temp_..C.=="26")
S2wRi26C <- subset(data, cell_line=="S2wRi" & temp_..C.=="26")

Dsim6B26C <- subset(data, cell_line=="Dsim6B" & temp_..C.=="26")
Dsim6BwMel26C <- subset(data, cell_line=="Dsim6BwMel" & temp_..C.=="26")
Dsim6BwWil26C <- subset(data, cell_line=="Dsim6BwWil" & temp_..C.=="26")


###############
### average and sd doubling times
##############
###################
###############
### average and sd doubling times
##############

mean(JW18DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(JW18DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(JW18wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(JW18wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(JW18wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(JW18wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(JW18wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(JW18wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(S2DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(S2DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(S2wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(S2wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(S2wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(S2wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(S2wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(S2wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(Dsim6B23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(Dsim6B23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(Dsim6BwMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(Dsim6BwMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

mean(Dsim6BwWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
sd(Dsim6BwWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

####################
##### average and sd growth rate per day
mean(JW18DOX23C$growth_rate_per_day)
sd(JW18DOX23C$growth_rate_per_day)

mean(JW18wMel23C$growth_rate_per_day)
sd(JW18wMel23C$growth_rate_per_day)

mean(JW18wRi23C$growth_rate_per_day)
sd(JW18wRi23C$growth_rate_per_day)

mean(JW18wWil23C$growth_rate_per_day)
sd(JW18wWil23C$growth_rate_per_day)

mean(S2DOX23C$growth_rate_per_day)
sd(S2DOX23C$growth_rate_per_day)

mean(S2wMel23C$growth_rate_per_day)
sd(S2wMel23C$growth_rate_per_day)

mean(S2wRi23C$growth_rate_per_day)
sd(S2wRi23C$growth_rate_per_day)

mean(S2wWil23C$growth_rate_per_day)
sd(S2wWil23C$growth_rate_per_day)

mean(Dsim6B23C$growth_rate_per_day)
sd(Dsim6B23C$growth_rate_per_day)

mean(Dsim6BwMel23C$growth_rate_per_day)
sd(Dsim6BwMel23C$growth_rate_per_day)

mean(Dsim6BwWil23C$growth_rate_per_day)
sd(Dsim6BwWil23C$growth_rate_per_day)

#####################
##### Infected vs uninfected
#####################
wilcox.test(JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(JW18DOX$growth_rate_per_day, JW18wMel$growth_rate_per_day)
wilcox.test(JW18DOX$cells.day, JW18wMel$cells.day)

> wilcox.test(JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and JW18wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 0, p-value = 0.0003108
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(JW18DOX$growth_rate_per_day, JW18wMel$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  JW18DOX$growth_rate_per_day and JW18wMel$growth_rate_per_day
W = 56, p-value = 0.001446
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(JW18DOX$growth_rate_per_day, JW18wMel$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(JW18DOX$cells.day, JW18wMel$cells.day)

	Wilcoxon rank sum test with continuity correction

data:  JW18DOX$cells.day and JW18wMel$cells.day
W = 56, p-value = 0.001446
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(JW18DOX$cells.day, JW18wMel$cells.day) :
  cannot compute exact p-value with ties
>


############
wilcox.test(JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(JW18DOX$growth_rate_per_day, JW18wRi$growth_rate_per_day)
wilcox.test(JW18DOX$cells.day, JW18wRi$cells.day)

> wilcox.test(JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and JW18wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 0, p-value = 0.006061
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(JW18DOX$growth_rate_per_day, JW18wRi$growth_rate_per_day)

	Wilcoxon rank sum exact test

data:  JW18DOX$growth_rate_per_day and JW18wRi$growth_rate_per_day
W = 28, p-value = 0.006061
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(JW18DOX$cells.day, JW18wRi$cells.day)

	Wilcoxon rank sum test with continuity correction

data:  JW18DOX$cells.day and JW18wRi$cells.day
W = 28, p-value = 0.01056
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(JW18DOX$cells.day, JW18wRi$cells.day) :
  cannot compute exact p-value with ties

###########

wilcox.test(JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(JW18DOX$growth_rate_per_day, JW18wWil$growth_rate_per_day)
wilcox.test(JW18DOX$cells.day, JW18wWil$cells.day)

> wilcox.test(JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  JW18DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and JW18wWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 4, p-value = 0.0002598
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(JW18DOX$growth_rate_per_day, JW18wWil$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  JW18DOX$growth_rate_per_day and JW18wWil$growth_rate_per_day
W = 86, p-value = 0.0009339
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(JW18DOX$growth_rate_per_day, JW18wWil$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(JW18DOX$cells.day, JW18wWil$cells.day)

	Wilcoxon rank sum test with continuity correction

data:  JW18DOX$cells.day and JW18wWil$cells.day
W = 87, p-value = 0.0006991
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(JW18DOX$cells.day, JW18wWil$cells.day) :
  cannot compute exact p-value with ties

##########
wilcox.test(S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(S2DOX$growth_rate_per_day, S2wMel$growth_rate_per_day)
wilcox.test(S2DOX$cells.day, S2wMel$cells.day)

> wilcox.test(S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and S2wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 0, p-value = 0.0001554
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(S2DOX$growth_rate_per_day, S2wMel$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  S2DOX$growth_rate_per_day and S2wMel$growth_rate_per_day
W = 64, p-value = 0.000931
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(S2DOX$growth_rate_per_day, S2wMel$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(S2DOX$cells.day, S2wMel$cells.day)

	Wilcoxon rank sum exact test

data:  S2DOX$cells.day and S2wMel$cells.day
W = 58, p-value = 0.004662
alternative hypothesis: true location shift is not equal to 0


##########
wilcox.test(S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(S2DOX$growth_rate_per_day, S2wRi$growth_rate_per_day)
wilcox.test(S2DOX$cells.day, S2wRi$cells.day)

> wilcox.test(S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and S2wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 1, p-value = 0.0003108
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(S2DOX$growth_rate_per_day, S2wRi$growth_rate_per_day)

	Wilcoxon rank sum exact test

data:  S2DOX$growth_rate_per_day and S2wRi$growth_rate_per_day
W = 63, p-value = 0.0003108
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(S2DOX$cells.day, S2wRi$cells.day)

	Wilcoxon rank sum exact test

data:  S2DOX$cells.day and S2wRi$cells.day
W = 42, p-value = 0.3282
alternative hypothesis: true location shift is not equal to 0

#####

##########
wilcox.test(S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(S2DOX$growth_rate_per_day, S2wWil$growth_rate_per_day)
wilcox.test(S2DOX$cells.day, S2wWil$cells.day)

wilcox.test(S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  S2DOX$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and S2wWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 0, p-value = 1.083e-05
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(S2DOX$growth_rate_per_day, S2wWil$growth_rate_per_day)

	Wilcoxon rank sum exact test

data:  S2DOX$growth_rate_per_day and S2wWil$growth_rate_per_day
W = 100, p-value = 1.083e-05
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(S2DOX$cells.day, S2wWil$cells.day)

	Wilcoxon rank sum exact test

data:  S2DOX$cells.day and S2wWil$cells.day
W = 69, p-value = 0.1655
alternative hypothesis: true location shift is not equal to 0

##########
wilcox.test(Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., Dsim6BwMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(Dsim6B$growth_rate_per_day, Dsim6BwMel$growth_rate_per_day)
wilcox.test(Dsim6B$cells.day, Dsim6BwMel$cells.day)

> wilcox.test(Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., Dsim6BwMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and Dsim6BwMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 25, p-value = 0.5054
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(Dsim6B$growth_rate_per_day, Dsim6BwMel$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  Dsim6B$growth_rate_per_day and Dsim6BwMel$growth_rate_per_day
W = 39, p-value = 0.4945
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(Dsim6B$growth_rate_per_day, Dsim6BwMel$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(Dsim6B$cells.day, Dsim6BwMel$cells.day)

	Wilcoxon rank sum exact test

data:  Dsim6B$cells.day and Dsim6BwMel$cells.day
W = 62, p-value = 0.0006216
alternative hypothesis: true location shift is not equal to 0



################################
wilcox.test(Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., Dsim6BwWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(Dsim6B$growth_rate_per_day, Dsim6BwWil$growth_rate_per_day)
wilcox.test(Dsim6B$cells.day, Dsim6BwWil$cells.day)

> wilcox.test(Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., Dsim6BwWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum test with continuity correction

data:  Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and Dsim6BwWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 120, p-value = 0.3566
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(Dsim6B$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.,  :
  cannot compute exact p-value with ties
> wilcox.test(Dsim6B$growth_rate_per_day, Dsim6BwWil$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  Dsim6B$growth_rate_per_day and Dsim6BwWil$growth_rate_per_day
W = 179, p-value = 0.3728
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(Dsim6B$growth_rate_per_day, Dsim6BwWil$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(Dsim6B$cells.day, Dsim6BwWil$cells.day)

	Wilcoxon rank sum test with continuity correction

data:  Dsim6B$cells.day and Dsim6BwWil$cells.day
W = 185, p-value = 0.281
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(Dsim6B$cells.day, Dsim6BwWil$cells.day) :
  cannot compute exact p-value with ties

##########
## wMel vs wRi
##########

wilcox.test(JW18wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(JW18wRi$growth_rate_per_day, JW18wMel$growth_rate_per_day)
wilcox.test(JW18wRi$cells.day, JW18wMel$cells.day)

> wilcox.test(JW18wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., JW18wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  JW18wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and JW18wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 32, p-value = 0.00404
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(JW18wRi$growth_rate_per_day, JW18wMel$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  JW18wRi$growth_rate_per_day and JW18wMel$growth_rate_per_day
W = 0, p-value = 0.00836
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(JW18wRi$growth_rate_per_day, JW18wMel$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(JW18wRi$cells.day, JW18wMel$cells.day)

	Wilcoxon rank sum exact test

data:  JW18wRi$cells.day and JW18wMel$cells.day
W = 0, p-value = 0.00404
alternative hypothesis: true location shift is not equal to 0

#########

wilcox.test(S2wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(S2wRi$growth_rate_per_day, S2wMel$growth_rate_per_day)
wilcox.test(S2wRi$cells.day, S2wMel$cells.day)

## ns ##

> wilcox.test(S2wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., S2wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  S2wRi$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and S2wMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 23, p-value = 0.3823
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(S2wRi$growth_rate_per_day, S2wMel$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  S2wRi$growth_rate_per_day and S2wMel$growth_rate_per_day
W = 41, p-value = 0.3717
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(S2wRi$growth_rate_per_day, S2wMel$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(S2wRi$cells.day, S2wMel$cells.day)

	Wilcoxon rank sum exact test

data:  S2wRi$cells.day and S2wMel$cells.day
W = 52, p-value = 0.03792
alternative hypothesis: true location shift is not equal to 0

#########

wilcox.test(Dsim6BwWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., Dsim6BwMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
wilcox.test(Dsim6BwWil$growth_rate_per_day, Dsim6BwMel$growth_rate_per_day)
wilcox.test(Dsim6BwWil$cells.day, Dsim6BwMel$cells.day)

## ns ##

> wilcox.test(Dsim6BwWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1., Dsim6BwMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)

	Wilcoxon rank sum exact test

data:  Dsim6BwWil$doubling_time_Td..t2.t1..ln.2..ln.q2.q1. and Dsim6BwMel$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.
W = 37, p-value = 0.6454
alternative hypothesis: true location shift is not equal to 0

> wilcox.test(Dsim6BwWil$growth_rate_per_day, Dsim6BwMel$growth_rate_per_day)

	Wilcoxon rank sum test with continuity correction

data:  Dsim6BwWil$growth_rate_per_day and Dsim6BwMel$growth_rate_per_day
W = 27.5, p-value = 0.674
alternative hypothesis: true location shift is not equal to 0

Warning message:
In wilcox.test.default(Dsim6BwWil$growth_rate_per_day, Dsim6BwMel$growth_rate_per_day) :
  cannot compute exact p-value with ties
> wilcox.test(Dsim6BwWil$cells.day, Dsim6BwMel$cells.day)

	Wilcoxon rank sum exact test

data:  Dsim6BwWil$cells.day and Dsim6BwMel$cells.day
W = 33, p-value = 0.9591
alternative hypothesis: true location shift is not equal to 0


##################################

> mean(JW18DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 2.219795
> sd(JW18DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.2742053
>
> mean(JW18wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 3.695692
> sd(JW18wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.3597444
>
> mean(JW18wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 21.03906
> sd(JW18wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 21.44571
>
> mean(JW18wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 4.149318
> sd(JW18wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.5934905

> mean(S2DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 2.130694
> sd(S2DOX23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.3936941
>
> mean(S2wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 4.344632
> sd(S2wMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.7820395
>
> mean(S2wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 3.589388
> sd(S2wRi23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.4518372
>
> mean(S2wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 3.335717
> sd(S2wWil23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.4065603

> mean(Dsim6B23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 3.325743
> sd(Dsim6B23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 0.7657272
>
> mean(Dsim6BwMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 4.846998
> sd(Dsim6BwMel23C$doubling_time_Td..t2.t1..ln.2..ln.q2.q1.)
[1] 2.09333
>
> ####################
> ##### average and sd growth rate per day
> mean(JW18DOX23C$growth_rate_per_day)
[1] 1.216
> sd(JW18DOX23C$growth_rate_per_day)
[1] 0.4457914
>
> mean(JW18wMel23C$growth_rate_per_day)
[1] 0.398
> sd(JW18wMel23C$growth_rate_per_day)
[1] 0.06760178
>
> mean(JW18wRi23C$growth_rate_per_day)
[1] 0.08333333
> sd(JW18wRi23C$growth_rate_per_day)
[1] 0.06027714
>
> mean(JW18wWil23C$growth_rate_per_day)
[1] 0.332
> sd(JW18wWil23C$growth_rate_per_day)
[1] 0.0858487

> mean(S2DOX23C$growth_rate_per_day)
[1] 1.568
> sd(S2DOX23C$growth_rate_per_day)
[1] 1.123063
>
> mean(S2wMel23C$growth_rate_per_day)
[1] 0.316
> sd(S2wMel23C$growth_rate_per_day)
[1] 0.1121606
>
> mean(S2wRi23C$growth_rate_per_day)
[1] 0.428
> sd(S2wRi23C$growth_rate_per_day)
[1] 0.1229634
>
> mean(S2wWil23C$growth_rate_per_day)
[1] 0.484
> sd(S2wWil23C$growth_rate_per_day)
[1] 0.09071935
>
> mean(Dsim6B23C$growth_rate_per_day)
[1] 0.564
> sd(Dsim6B23C$growth_rate_per_day)
[1] 0.3203592
>
> mean(Dsim6BwMel23C$growth_rate_per_day)
[1] 0.458
> sd(Dsim6BwMel23C$growth_rate_per_day)
[1] 0.5048465



> ###########################
> ### starting concentration
> ###########################
> mean(JW18DOX23C$starting_concentration_.cells.mL.)
[1] 904600
> sd(JW18DOX23C$starting_concentration_.cells.mL.)
[1] 246689.7
>
> mean(JW18wMel23C$starting_concentration_.cells.mL.)
[1] 1066600
> sd(JW18wMel23C$starting_concentration_.cells.mL.)
[1] 349267.5
>
> mean(JW18wRi23C$starting_concentration_.cells.mL.)
[1] 761000
> sd(JW18wRi23C$starting_concentration_.cells.mL.)
[1] 178530.1
>
> mean(JW18wWil23C$starting_concentration_.cells.mL.)
[1] 1480000
> sd(JW18wWil23C$starting_concentration_.cells.mL.)
[1] 497644.5

> mean(S2DOX23C$starting_concentration_.cells.mL.)
[1] 826600
> sd(S2DOX23C$starting_concentration_.cells.mL.)
[1] 350831.6
>
> mean(S2wMel23C$starting_concentration_.cells.mL.)
[1] 1500000
> sd(S2wMel23C$starting_concentration_.cells.mL.)
[1] 430871.2
>
> mean(S2wRi23C$starting_concentration_.cells.mL.)
[1] 2134000
> sd(S2wRi23C$starting_concentration_.cells.mL.)
[1] 980244.9
>
> mean(S2wWil23C$starting_concentration_.cells.mL.)
[1] 1866000
> sd(S2wWil23C$starting_concentration_.cells.mL.)
[1] 639515.4

> mean(Dsim6B23C$starting_concentration_.cells.mL.)
[1] 6052000
> sd(Dsim6B23C$starting_concentration_.cells.mL.)
[1] 3050798
>
> mean(Dsim6BwMel23C$starting_concentration_.cells.mL.)
[1] 2684000
> sd(Dsim6BwMel23C$starting_concentration_.cells.mL.)
[1] 1197531

> mean(Dsim6BwWil23C$starting_concentration_.cells.mL.)
[1] 2392000
> sd(Dsim6BwWil23C$starting_concentration_.cells.mL.)
[1] 1393223

> ##################################
> ### ending concentration
> ##################################
> mean(JW18DOX23C$ending_concentration_.cells.mL.)
[1] 8008000
> sd(JW18DOX23C$ending_concentration_.cells.mL.)
[1] 396383.7
>
> mean(JW18wMel23C$ending_concentration_.cells.mL.)
[1] 3928000
> sd(JW18wMel23C$ending_concentration_.cells.mL.)
[1] 925699.7
>
> mean(JW18wRi23C$ending_concentration_.cells.mL.)
[1] 1136667
> sd(JW18wRi23C$ending_concentration_.cells.mL.)
[1] 66583.28
>
> mean(JW18wWil23C$ending_concentration_.cells.mL.)
[1] 4824000
> sd(JW18wWil23C$ending_concentration_.cells.mL.)
[1] 1371598

> mean(S2DOX23C$ending_concentration_.cells.mL.)
[1] 7920000
> sd(S2DOX23C$ending_concentration_.cells.mL.)
[1] 1717265
>
> mean(S2wMel23C$ending_concentration_.cells.mL.)
[1] 4764000
> sd(S2wMel23C$ending_concentration_.cells.mL.)
[1] 1305328
>
> mean(S2wRi23C$ending_concentration_.cells.mL.)
[1] 8152000
> sd(S2wRi23C$ending_concentration_.cells.mL.)
[1] 3023660
>
> mean(S2wWil23C$ending_concentration_.cells.mL.)
[1] 8004000
> sd(S2wWil23C$ending_concentration_.cells.mL.)
[1] 2263245
>

> mean(Dsim6B23C$ending_concentration_.cells.mL.)
[1] 24600000
> sd(Dsim6B23C$ending_concentration_.cells.mL.)
[1] 9264718
>
> mean(Dsim6BwMel23C$ending_concentration_.cells.mL.)
[1] 8372000
> sd(Dsim6BwMel23C$ending_concentration_.cells.mL.)
[1] 1723171

>
> mean(Dsim6BwWil23C$ending_concentration_.cells.mL.)
[1] 8312000
> sd(Dsim6BwWil23C$ending_concentration_.cells.mL.)
[1] 2509526
>
