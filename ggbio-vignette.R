# GGBIO TUTORIALS
# 20150108
# 
library(ggbio)
p.ideo <- Ideogram(genome = "hg19")
p.ideo
library(GenomicRanges)
## special highlights instead of zoomin!
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))
# 
# library(OrganismDbi)
library(ggbio)
library(Homo.sapiens)
class(Homo.sapiens)
##
data(genesymbol, package = "biovizBase")
genesymbol

wh <- genesymbol[c("BRCA1", "NBR1")]
wh
wh <- range(wh, ignore.strand = TRUE)
p.txdb <- autoplot(Homo.sapiens, which = wh)
p.txdb
autoplot(Homo.sapiens, which = wh, label.color = "black", color = "brown", fill = "brown")

# To change the intron geometry, use gap.geom to control it, check out geom alignment for more control parameters.
autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")

# To collapse all features, use stat ’reduce’
autoplot(Homo.sapiens, which = wh, stat = "reduce")

# Label could be turned off by setting it to FALSE, you could also use expression to make a flexible label combination from column names.
columns(Homo.sapiens)
autoplot(Homo.sapiens, which = wh, columns = c("TXNAME", "GO"), names.expr = "TXNAME::GO")

# TxDb doesn’t contain any gene symbol information, so we use tx id as default for label
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
autoplot(txdb, which = wh)

# 
library(biovizBase)
gr.txdb <- crunch(txdb, which = wh)
?crunch
## change column to ✬model✬
colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$tx_id)
## fake some randome names
names(grl) <- sample(LETTERS, size = length(grl), replace = TRUE)
grl

# We get our example data ready, it meets all requirements, to make it a gene model track it’s pretty simple to use autoplot,
# but don’t forget mapping because we changed our column names, asssume you store you model key words in column
# ’model’.
autoplot(grl, aes(type = model))
ggplot() + geom_alignment(grl, type = "model")

# To add a reference track, we need to load a BSgenome object from the annotation package. You can choose to plot the sequence as text, rect, segment
# You can pass a zoom in factor into zoom function, if it’s over 1 it’s zooming out, if it’s smaller than 1 it’s zooming in.
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
p.bg <- autoplot(bg, which = wh)
## no geom
p.bg
## segment
p.bg + zoom(1/100)
## rectangle
p.bg + zoom(1/1000)
## text
p.bg + zoom(1/2500)
# To override a zemantic zoom threshold, you simply provide a geom explicitly.
library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
## force to use geom ✬segment✬ at this level
autoplot(bg, which = resize(wh, width = width(wh)/2000), geom = "segment")
