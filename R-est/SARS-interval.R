# Explore estimate of SARS serial interval, used by Kissler, et al., taken
# from Lipsitch, et al. "Transmission Dynamics and Control of Severe Acute
# Respiratory Syndrome".  Data from Fig. 1 (E) in that paper.

data <- c (3, 1, 14, 15, 10, 11, 19, 21, 22, 18, 9, 12, 8, 6, 3, 2, 3, 1, 1, 1)

SARS_shape <- 2.35
SARS_scale <- 9.48

SARS_gen_interval <- dweibull (0:21, shape=SARS_shape, scale=SARS_scale)

plot(data, pch=20, ylim=c(0,1.05*max(data)), yaxs="i", 
     xlim=c(0,length(data)+1), xaxs="i",
     xlab="", ylab="", type="n")

abline(h=c(0,5,10,15,20))

lines (0:(length(SARS_gen_interval)-1),
       SARS_gen_interval*sum(data),
       lwd=3, col="gray")

for (i in 1:length(data)) lines (c(i,i),c(0,data[i]),lwd=2)
