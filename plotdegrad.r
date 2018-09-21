#plotdegrad is a function used to plot area with rectangles of grey shades
#author: Diane Bailleul
#email: diane.bailleul.pro@gmail.com

plotdegrad <- function(nb_degrad = NULL, attenua = NULL, xleft0 = NULL, ybottom0 = NULL, xright0 = NULL, ytop0 = NULL){
	a <- nb_degrad
	dimenx <- xright0/a
	b <- a+attenua
	greycol <- grey(0:b/b)
	greycol <- greycol[-c(1:attenua)]
		for (i in 1:a){
			rect(xleft = xleft0 +(i-1)*dimenx, ybottom = ybottom0, xright = xright0-(a-i)*dimenx, ytop = ytop0,
			border = FALSE, col = greycol[i])
	}
}
