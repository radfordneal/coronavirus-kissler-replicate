# Explore estimate of SARS serial interval, used by Kissler, et al., taken
# from Lipsitch, et al. "Transmission Dynamics and Control of Severe Acute
# Respiratory Syndrome".  Data from Fig. 1 (E) in that paper.

pdf ("SARS-interval.pdf", height=8, width=6)
par(mar=c(1.5,2.3,3,0.5),mgp=c(1.4,0.3,0),tcl=-0.22)

data <- c (3, 1, 14, 15, 10, 11, 19, 21, 22, 18, 9, 12, 8, 6, 3, 2, 3, 1, 1, 1)

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <- dweibull (0:21, shape=SARS_shape, scale=SARS_scale)

par(mfrow=c(3,1))

plot(data, pch=20, ylim=c(0,1.05*max(data)), yaxs="i", 
     xlim=c(0,length(data)+1), xaxs="i",
     xlab="", ylab="", type="n")

abline(h=c(0,5,10,15,20))

lines (0:(length(SARS_gen_interval)-1),
       SARS_gen_interval*sum(data),
       lwd=3, col="gray")

for (i in 1:length(data)) lines (c(i,i),c(0,data[i]),lwd=2)

title (
 "Data on SARS serial interval distribution and Weibull fit, from Lipsitch, et al.")

dev.off()
