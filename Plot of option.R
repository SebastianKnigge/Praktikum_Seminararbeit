###################################
#  plot der put option Auszahlung #

x <- seq(0,2, length.out = 100000)
K <- 1.1

y <- ifelse(K-x>0,K-x,0)

plot(x, y, type="l", ylim=c(0,1.1), lwd=2, ylab="option price", xlab="underlying")
abline(v=1.1, col=2)
legend(x="topright", legend=c("price of the option", "strike price = 1.1"), col=1:2, lwd=1, cex=0.5, inset=0.01)
