
#Read data
#Data is copy pasted from matlab struct for test
eog_grads <- read.csv("eog_grads.csv")

#Exclude outlier (#4)
eog.ex = eog_grads$eog[c(1:3, 5:22)]
grads.ex = eog_grads$grads[c(1:3, 5:22)]

svg("../../Analysis Output/eog_grads.svg")

plot(eog_grads$eog, eog_grads$grads,
  col = ifelse(1:nrow(eog_grads) == 4, "red", "black"),
  
main = 'Gradiometer amp ~ EOG amp',
ylab = 'Gradiometer amplitude',
xlab = 'EOG amplitude')

abline(lm(eog_grads$grads ~ eog_grads$eog), col = "red")
abline(lm(grads.ex ~ eog.ex))

text(-0.00025, 1.36e-11, expression(paste(R^2==-0.44, ",  p = 0.04")), col = 'red', pos = 4)
text(-0.00025, 1.3e-11, expression(paste(R^2==0.06, ",  p = 0.79")), col = 'black', pos = 4)

cor.test(eog_grads$grads, eog_grads$eog, method = "pearson")
cor.test(grads.ex, eog.ex, method = "pearson")

dev.off()