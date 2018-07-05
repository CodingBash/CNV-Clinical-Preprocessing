install.packages("CNprep")
library(CNprep)
CNprep::

data(segexample)
data(ratexample)
data(normsegs)
#small toy example
segtable<-CNpreprocessing(segall=segexample[segexample[,"ID"]=="WZ1",],
                          ratall=ratexample,"ID","start","end",chromcol="chrom",bpstartcol="chrom.pos.start",
                          bpendcol="chrom.pos.end",blsize=50,minjoin=0.25,cweight=0.4,bstimes=50,
                          chromrange=1:3,distrib="Rparallel",njobs=2,modelNames="E",
                          normalength=normsegs[,1],normalmedian=normsegs[,2])
## Not run:
#Example 1: 5 whole genome analysis, choosing the right format of arguments
segtable<-CNpreprocessing(segall=segexample,ratall=ratexample,"ID","start","end",
                          chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=50,
                          minjoin=0.25,cweight=0.4,bstimes=50,chromrange=1:22,distrib="Rparallel",njobs=40,
                          modelNames="E",normalength=normsegs[,1],normalmedian=normsegs[,2])
#Example 2: how to use annotexample, when segment table does not have columns of
#integer postions in terms of measuring units(probes), such as "mysegs" below
mysegs<-segexample[,c(1,5:12)]
data(annotexample)
segtable<-CNpreprocessing(segall=mysegs,ratall=ratexample,"ID",chromcol="chrom",
                          bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",annot=annotexample,
                          annotstartcol="CHROM.POS",annotendcol="CHROM.POS",annotchromcol="CHROM",
                          blsize=50,minjoin=0.25,cweight=0.4,bstimes=50,chromrange=1:22,distrib="Rparallel",
                          njobs=40,modelNames="E",normalength=normsegs[,1],normalmedian=normsegs[,2])
## End(Not run)
