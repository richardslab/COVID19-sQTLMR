library(data.table)
library(grid)
library(gridExtra)
library(lattice)
library(plotrix)

library(plotrix)



DT_pp <- fread("pptable.csv", header = TRUE)
colnames(DT_pp) <- c("", "A2", "B2", "C2")

DT_A2 <- fread("A2_OAS1_p_snv.txt")
DT_A2[, mlog10p.o := -log10(P_A2)]
DT_A2[, mlog10p.e := -log10(P_OAS1)]
DT_A2[, mbp := POS/1000000]
DT_A2[, r2bin := .bincode(R2, seq(0,1,by=0.2), right = TRUE, include.lowest = TRUE)]

DT_B2 <- fread("B2_eur_OAS1_p_snv.txt")
DT_B2[, mlog10p.o := -log10(P_B2_eur)]
DT_B2[, mlog10p.e := -log10(P_OAS1)]
DT_B2[, mbp := POS/1000000]
DT_B2[, r2bin := .bincode(R2, seq(0,1,by=0.2), right = TRUE, include.lowest = TRUE)]

DT_C2 <- fread("C2_eur_OAS1_p_snv.txt")
DT_C2[, mlog10p.o := -log10(P_C2_eur)]
DT_C2[, mlog10p.e := -log10(P_OAS1)]
DT_C2[, mbp := POS/1000000]
DT_C2[, r2bin := .bincode(R2, seq(0,1,by=0.2), right = TRUE, include.lowest = TRUE)]

png(file="oas1-covid19-coloc.png", width=1800*2, height=2700*2, res=600)

# colbin <- colorRampPalette(c("#d0d1e6", "#023858"))(10)
colbin <- rev(c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#357EBD"))


par(mar=c(2,4,1,1))

layout(mat = matrix(c(1:10), nrow = 5, ncol = 2),
        heights=c(0.15,0.5,0.5,0.5,0.5,
                  0.15,0.5,0.5,0.5,0.5),
        widths=c(0.5,0.5,0.5,0.5,0.5,
                 0.5,0.5,0.5,0.5,0.5)
       #       heights = c(1,1,1,1,1), # Heights of the two rows
#       widths = c(1,1)
) # Widths of the two columns

par(mar=c(2,4,1,1))
plot.new()
DT_A2[, plot(mlog10p.e,mlog10p.o, xlab="", xaxt="n", ylab=expression(A2~~-log[10](italic(p))), bty="n", bg=colbin[r2bin], pch=21, lwd=0.7)]
DT_A2[SNP=="rs4767027", points(mlog10p.e,mlog10p.o, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]

DT_B2[, plot(mlog10p.e,mlog10p.o, xlab="", xaxt="n", ylab=expression(B2~~-log[10](italic(p))), bty="n", bg=colbin[r2bin], pch=21, lwd=0.7)]
DT_B2[SNP=="rs4767027", points(mlog10p.e,mlog10p.o, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]

DT_C2[, plot(mlog10p.e,mlog10p.o, xlab=expression(OAS1~~-log[10](italic(p))), ylab=expression(C2~~-log[10](italic(p))), bty="n", bg=colbin[r2bin], pch=21, lwd=0.7)]
DT_C2[SNP=="rs4767027", points(mlog10p.e,mlog10p.o, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]
mtext(expression(OAS1~~-log[10](italic(p))), side=1, line = 3, cex = 0.7)

plot(1,1, bty="n", col="white", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,2), ylim=c(0,2))
addtable2plot(x=0.13, y=0.75,table=DT_pp, xjust=0, yjust=0.5, cex=1.1, hlines = TRUE, vlines = TRUE, xpad=0.4, ypad=0.8, title="Coloc Posterior Probability")

par(mar=c(1,15,3,1))
lut=colbin
min=0
max=1
scale = (length(lut))/(max-min)
nticks=6
ticks=seq(min, max, len=nticks)
title=''
plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
mtext(expression(paste( italic(r^2), " to lead SNP")), side=2, line=0, cex=0.5, las=1)
axis(3, ticks, las=1, line = 0, mgp=c(0,0.7,0), cex.axis=0.6)
for (i in 1:(length(lut))) {
        x = (i-1)/scale + min
        rect(x,0,x+1/scale, 10,col=lut[i], border=1)
}

par(mar=c(2,3,1,1))
DT_A2[, plot(mbp,mlog10p.o, xlab="", xaxt="n", ylab="", bty="n", bg=colbin[r2bin], pch=21, lwd=0.7)]
DT_A2[SNP=="rs4767027", points(mbp,mlog10p.o, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]

#legend("topright",
#       legend=c(".8-1", ".6-.8", ".4-.6", ".2-.4", "0-.2"),
#       fill=rev(colbin),
#       cex=0.7,
#       y.intersp=0.8,
#       x.intersp=0.2,
#       inset=0.01,
#       border = rev(colbin),
#       box.col = "white",
#       title = expression( italic(r^2)))


#lgd_ = rep(NA, 5)
#lgd_[c(1,2,3,4,5)] = c("0.8-1","0.6-0.8","0.4-0.6","0.2-0.4", "0-0.2")
#legend(x = 113.55, y = 9.5, bty = "n",
#       legend = lgd_,
#       fill = rev(colbin),
#       border = "black",
#       y.intersp = 1,
#       x.intersp = 0.25,
#       cex = 1)
#text(113.67,9.6,labels = expression('R'^"2"), cex=1.5)

legend_image <- as.raster(matrix(colbin, ncol=1))

DT_B2[, plot(mbp,mlog10p.o, xlab="", xaxt="n", ylab="", bty="n", bg=colbin[r2bin], pch=21, lwd=0.7)]
DT_B2[SNP=="rs4767027", points(mbp,mlog10p.o, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]
DT_C2[, plot(mbp,mlog10p.o, xlab="", xaxt="n", ylab="", bty="n", bg=colbin[r2bin], pch=21, lwd=0.7)]
DT_C2[SNP=="rs4767027", points(mbp,mlog10p.o, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]

par(mar=c(4,3,1,1))
DT_C2[, plot(mbp,mlog10p.e, xlab="Chromosome 12 (Mbp)",  xaxt="n", bty="n", bg=colbin[r2bin], pch=21, lwd=0.7, ylab="")]
DT_C2[SNP=="rs4767027", points(mbp,mlog10p.e, pch=23, bg="#9632B8", lwd=0.7, cex=1.3)]
axis(1,at=round(seq(min(DT_C2$mbp), max(DT_C2$mbp), by=0.2),digits = 1), labels = round(seq(min(DT_C2$mbp), max(DT_C2$mbp), by=0.2),digits = 1))
mtext(expression(OAS1~~-log[10](italic(p))), side=2, line = 3, cex = 0.7)
text(x=DT_C2[SNP=="rs4767027", mbp], y=DT_C2[SNP=="rs4767027", mlog10p.e], labels = "rs4767027", col="black", pos = 4, offset=0.5)
dev.off()

