args <- commandArgs(TRUE)
# manhattan (raleigh) plot function
manhattan.plot <- function(pval = pval, snp = snp, low.col = "lightblue", high.col = "steelblue", threshold = 1e-5, ...) {
  # set chromosome length and gaps
  chr.len <- c(22422827, 23011544, 21146708, 24543557, 27905053, 1351857)
  names(chr.len) <- c("X", "2L", "2R", "3L", "3R", "4")
  chr.gap <- 2000000
  chr.cum <- cumsum(c(0, chr.len[1:5])) + (0:5)*chr.gap
  names(chr.cum) <- names(chr.len)
  # extract necessary chromosome coordinates
  chr.info <- unlist(strsplit(snp, split = "_"))
  chr <- chr.info[seq(1, length(chr.info), 3)]
  pos <- as.numeric(chr.info[seq(2, length(chr.info), 3)])
  # pval (take -log10)
  pval <- -log10(pval)
  
  # get rid of chr 4 if it's not in the data set
  if (sum(chr == "4") == 0) {
    chr.len <- chr.len[c("X", "2L", "2R", "3L", "3R")]
    chr.cum <- chr.cum[c("X", "2L", "2R", "3L", "3R")]
  }
  # plot
  plot(chr.cum[chr] + pos, pval, col = ifelse(pval > -log10(threshold), high.col, low.col), pch = 16, axes = F, ylab = expression(paste("-log", {}[10], "(", italic(P)," value)")), xlab = "", ..., mgp = c(1.5, 0.5, 0), cex = 0.6)
  axis(side = 1, at = chr.cum + chr.len/2, labels = parse(text = paste("italic(\"", names(chr.len), "\")", sep = "")), mgp = c(2.5, 0.5, 0), tck = -0.015, cex.axis = 6/par("ps")/par("cex"))
  axis(side = 2, las = 1, mgp = c(1.5, 0.5, 0), tck = -0.015, cex.axis = 6/par("ps")/par("cex"))
  abline(h = -log10(threshold), lty = 2)
  box()
}

pval <- read.table(paste(args[1], "/gwas.all.assoc", sep = ""), header = TRUE, as.is = TRUE)
png(file = args[2], res = 300, width = 8, height = 8, unit = "in")
par(mfrow = c(4, 1))
manhattan.plot(pval = pval$FemaleMixedPval, snp = pval$ID, ylim = c(0, 10), main = paste(args[3], " Female", sep = ""))
manhattan.plot(pval = pval$MaleMixedPval, snp = pval$ID, ylim = c(0, 10), main = paste(args[3], " Male", sep = ""))
manhattan.plot(pval = pval$AvgMixedPval, snp = pval$ID, ylim = c(0, 10), main = paste(args[3], " Sex Average", sep = ""))
manhattan.plot(pval = pval$DiffMixedPval, snp = pval$ID, ylim = c(0, 10), main = paste(args[3], " Sex Difference", sep = ""))
dev.off()
