#This is a exemple of the R scripts used in the article "Large‐scale genetic 
##panmixia in the blue shark (Prionace glauca): A single worldwide population, 
##or a genetic lag‐time effect of the “grey zone” of differentiation?" to
##determine the existence and duration of the grey zone, for the graph generation.

#written by Diane Bailleul
#email: diane.bailleul.pro@gmail.com

#Please, do not use (with or witout modifications) without citing
##R and the original article.

##################################################
#data treatment and graph generation for N=100000#
##################################################

#rename data from SimuPop/grey_zone.py into res_simupopN100000.txt
data1 <- read.table("res_simupopN100000.txt", header = FALSE, sep="\t")

names(data1) <- c("gen", "Fst_tot", "rep", "Fst_sub", paste("sim", 1:1000, sep="_"))
head(data1)
data1 <- data1[,-1005]
data1$gen <- data1$gen + 20 #I don't like generation "0"

#pvalues computation
datat <- data1
pval_two <- NULL

for(i in 1:nrow(datat)){
	sim <- datat[i,5:1004]
	#pval_upper <- c(pval_upper, mean(sim <= datat[i,4]))
	pval_two <- c(pval_two, if (2*min(mean(sim <= datat[i,4]), mean(sim >= datat[i,4])) > 1) {1} else {2*min(mean(sim <= datat[i,4]), mean(sim >= datat[i,4]))})
}

#pval_upper <- 1 - pval_upper

#detection % per generation
ndata1 <- cbind(data1[,1:4], pval_two)
reccords <- unlist(lapply(split(ndata1[,5], ndata1[,1]), function(x) mean(x <= 0.05)))

Fst_tot <- unlist(lapply(split(ndata1[,2], ndata1[,1]),unique))

Fst_sub <- split(ndata1[,4], ndata1[,1])

#Fst_subb <- lapply(Fst_sub, function(x) x*10) #if necessary to change scale

#graph, at least
x <- as.numeric(names(reccords))
y <- as.vector(reccords)
lis <- loess(y ~ x, span = 0.5)
z   <- predict(lis, x)
preN100000 <- z
atN100000 <- x

par(mar = c(5.1, 4.1, 4.1, 4.1))

plot(Fst_tot~x, type = "l", ylim = c(0,1), xaxt = "n", yaxt = "n", 
	xlab = "", ylab = "", axes = FALSE)

mtext("Fst", side = 4, line = 3, par(las = 0), col = "darkblue")
mtext("generation", side = 1, line = 3, par(las = 0), col = "black")
mtext("detection capacity", side = 2, line = 3, par(las = 0), col = "darkgreen")

plotdegrad(nb_degrad = 2000, attenua = 700, xleft0 = 0, ybottom0 = 0, xright0 = 1170, ytop0 = 1)

#a <- 2000 #number of shades
#xright0 <- 250 #total length of the square
#xleft0 <- 0
#ybottom <- 0
#ytop0 <- 1

ytick <- seq(0, 1, by = 0.2)
axis(side = 2, at = ytick, labels = FALSE, col = "darkgreen")
text(x = -180, y = ytick + 0.03, col = "darkgreen", cex = 1,
     labels = c("0%","20%","40%","60%","80%","100%"), srt = 90, pos = 1, xpd = TRUE)

ztick <- seq(0, 1, by = 0.2)
axis(side = 4, at = ztick, labels = FALSE, col ="darkblue")
text(x = 2180, y = ztick + 0.02, col ="darkblue", cex = 1,
     labels = c("0","0.02","0.04","0.06","0.08","0.1"), srt = 270, pos = 1, xpd = TRUE)

xtick <- seq(0, 2000, by = 200)
axis(side = 1, at = xtick, labels = FALSE)
text(x = xtick, y = -0.1, 
     labels = xtick, srt = 0, xpd = TRUE, cex = 1)

lines(Fst_tot~x, type = "l", col = "darkblue", lwd = 2)
lines(as.vector(unlist(lapply(Fst_sub, median)))~x, type = "l", col = "royalblue")

lines(as.vector(unlist(lapply(Fst_sub, function(x) quantile(x, probs = 0.025))))~x, 
	lty = 3, col = "royalblue")
lines(as.vector(unlist(lapply(Fst_sub, function(x) quantile(x, probs = 0.975))))~x, 
	lty = 3, col = "royalblue")
lines(z~x, type = "l", col = "darkgreen", lwd = 2)

#Legends
#library(gplots)
legend(x = 280, y = 0.95,
            c("significant sub-Fst proportion", "Fst", "median sub-Fst", "sub-Fst 95% CI"),
            col = c("darkgreen", "darkblue", "royalblue", "royalblue"), 
		bg = c("white"), 
		#pch = c(16, 16),
		lwd = c(3, 3, 2, 2),
		lty = c(1, 1, 1, 3),
		inset = 0, 
		cex = 1,
		border = "white")

dev.copy2eps(file="RplotN100000.eps")
