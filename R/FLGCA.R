
fragmentListRead = function(file,nrows=-1) {
  read.table(file,header=T,nrows=nrows,comment.char="#");
}

fragmentListAnalyse = function(sampleName,fragmentListFile,qlim=c(0.25,0.75),downSample=NULL) {
  
  if (inherits(class(fragmentListFile),what="data.frame")) {
    flAll = fragmentListFile;
  } else {
    cat("Reading fragments file...",file=stderr()); flush(stderr());
    flAll = fragmentListRead(fragmentListFile);
    cat("done\n",file=stderr()); flush(stderr());
  }
  
  nnnReadFragments = flAll$N.LAST<=35 | flAll$N.FIRST <= 35;
  flAll = subset(flAll,subset=!nnnReadFragments);

  N.MAPPED = integer(nrow(flAll));
  GC.MAPPED = integer(nrow(flAll));
  N.MAPPED[!is.na(flAll$MQ.FIRST)] = flAll$N.FIRST[!is.na(flAll$MQ.FIRST)];
  GC.MAPPED[!is.na(flAll$MQ.FIRST)] = flAll$GC.FIRST[!is.na(flAll$MQ.FIRST)];
  N.MAPPED[!is.na(flAll$MQ.LAST)] = N.MAPPED[!is.na(flAll$MQ.LAST)] +  flAll$N.LAST[!is.na(flAll$MQ.LAST)];
  GC.MAPPED[!is.na(flAll$MQ.LAST)] = GC.MAPPED[!is.na(flAll$MQ.LAST)] + flAll$GC.LAST[!is.na(flAll$MQ.LAST)];
  flAll$N.MAPPED = N.MAPPED;
  flAll$GC.MAPPED = GC.MAPPED;
  fl = subset(flAll,subset= FLAGS == 0);
  count = nrow(fl);
  #if (!is.null(downSample)) {
  #  fl = fl[sample(count,downSample),];
  #}
  quartiles = quantile(fl$LENGTH,probs=c(0,0.25,0.50,0.75,1));
  median = quartiles[3];
  iqr = quartiles[4] - quartiles[3];
  lIqr= quartiles[4] - quartiles[2];
  lowEnd = median - lIqr * 2;
  highEnd = median + lIqr * 2;
  
  qlim.range = quantile(fl$LENGTH,probs=qlim);
  qlim.fl = fl[fl$LENGTH >= qlim.range[1] & fl$LENGTH >= qlim.range[2],];
  qlim.gc.pc = qlim.fl$GC.COUNT / qlim.fl$N.COUNT;
  cat("Fragments size vs gc correlation test...",file=stderr()); flush(stderr());  
  qlim.gc.cor.test = cor.test(qlim.gc.pc,qlim.fl$LENGTH);
  cat("done\n",file=stderr()); flush(stderr());
  cat("Fragments size mixture model fit... ",file=stderr()); flush(stderr());
  
  if (is.null(downSample)) downSample = min(nrow(fl),10000); 
  lengthsSample = if (nrow(fl) <= downSample) fl$LENGTH else sample(fl$LENGTH,downSample,replace=T);
    
  gaussFit = fragmentListFitGaussian(lengthsSample);
  cat("done\n",file=stderr()); flush(stderr());
  list(data.all=flAll,data=fl,quartiles=quartiles,median=median,mean=mean(qlim.fl$LENGTH,na.rm=T),sd=sd(qlim.fl$LENGTH,na.rm=T),
       qlim=qlim,qlim.range=qlim.range,qlim.gc.cor.test=qlim.gc.cor.test,count=count,gaussFit=gaussFit,
       n35.reads.count=length(which(nnnReadFragments)),lengthFitSample=lengthsSample,iqr=iqr,lowEnd=lowEnd,highEnd=highEnd);
}

fragmentListFitGaussian = function(lengths) {
   Mclust(data=lengths,modelNames="V");
}

regionBiasCounts = function(reference,seqName,startPos,endPos) {
  gc = length(which(reference.GC[[seqName]][startPos:endPos]));
  nuc = length(which(!referenceMask.N[[seqName]][startPos:endPos]));
  c(gc,nuc);
}

writeFragmentListResults = function(analysis,outDir,hist.binwidth=2,hist.sigma=4, hist.xlim=NULL,sampleName='') {
   histo.file = paste(outDir,"lengths.svg",sep="/");
   gc.file = paste(outDir,"gc.svg",sep="/");
   gc.length.file =paste(outDir,"gc-length.svg",sep="/");
   data.file = paste(outDir,"data.RData",sep="/");
   length.fit.plot.file = paste(outDir,"lengths-fit.pdf",sep="/");
   summary.file = paste(outDir,"summary.tsv",sep="/");
   if (is.null(downSample)) downSample = nrow(analysis$data);
   plot.data = if (nrow(analysis$data) <= downSample) analysis$data else analysis$data[sample.int(nrow(analysis$data),downSample),];
   
   cat("Writing statistics and graphics...",file=stderr()); flush(stderr());
   if (is.null(hist.xlim))
     hist.xlim = c(max(0,min(analysis$quartile[1],analysis$median - analysis$sd * hist.sigma)),
            max(analysis$median + analysis$sd * hist.sigma,analysis$quartiles[4]));
   
   pdf(file=length.fit.plot.file,width=10,height=4);
     layout(matrix(c(1,2,3),nrow=1));  
     lengths = analysis$lengthFitSample;
     plot(analysis$gaussFit,what= "BIC",main=sampleName);
     plot(analysis$gaussFit,what= "uncertainty",xlim=c(0,median(lengths,na.rm=T) + 4 * sd(lengths,na.rm=T)));
     plot(analysis$gaussFit,what= "density",xlim=c(0,median(lengths,na.rm=T) + 4 * sd(lengths,na.rm=T)));
   dev.off();
   
   lIqr=analysis$quartiles[4] - analysis$quartiles[2];
   plot.data$lowEnd = rep(analysis$median - lIqr * 2,nrow(plot.data));
   plot.data$highEnd = rep(analysis$median + lIqr * 2,nrow(plot.data));
   ggplot(plot.data, aes(x=LENGTH)) + coord_cartesian(xlim=hist.xlim) +
     geom_density(aes(y=..count..),alpha=.2,fill="#FF6666") +
     geom_vline(aes(xintercept=mean(LENGTH,na.rm=T)),color="red",linetype="dashed",size=1) +
     geom_vline(aes(xintercept=median(LENGTH,na.rm=T)),color="blue",linetype="solid",size=1) +
     geom_vline(aes(xintercept=mean(lowEnd)),color="black",linetype="solid",size=0.25) +
     geom_vline(aes(xintercept=mean(highEnd)),color="black",linetype="solid",size=0.25) +
     geom_vline(aes(xintercept=quantile(LENGTH,probs=0.25,na.rm=T)),color="blue",linetype="dashed",size=0.75) +
     geom_vline(aes(xintercept=quantile(LENGTH,probs=0.75,na.rm=T)),color="blue",linetype="dashed",size=0.75) +
     xlab("Fragment Length (bp)") + ylab("Frequency") +
     opts(title=sampleName);
   ggsave(histo.file,width=7,height=6,units="in");
   ggplot(plot.data, aes(x=GC.COUNT/N.COUNT)) + geom_density(aes(y=..count..),alpha=.2,fill="#FF6666") +
     geom_vline(aes(xintercept=mean(GC.COUNT/N.COUNT,na.rm=T)),color="red",linetype="dashed",size=1) +
     geom_vline(aes(xintercept=median(GC.COUNT/N.COUNT,na.rm=T)),color="blue",linetype="solid",size=1) +
     geom_vline(aes(xintercept=quantile(GC.COUNT/N.COUNT,probs=0.25,na.rm=T)),color="blue",linetype="dashed",size=0.75) +
     geom_vline(aes(xintercept=quantile(GC.COUNT/N.COUNT,probs=0.75,na.rm=T)),color="blue",linetype="dashed",size=0.75) +
     xlab("Fragment GC content (fraction)") + ylab("Frequency") +
     opts(title=sampleName);
   ggsave(gc.file,width=7,height=6,units="in");
   ggplot(plot.data) + geom_density2d(aes(x=LENGTH,y=GC.COUNT/N.COUNT),color="black",fill="#FF6666") +
     geom_vline(aes(xintercept=median(LENGTH,na.rm=T)),color="blue",linetype="solid",size=1) +
     geom_hline(aes(yintercept=median(GC.COUNT/N.COUNT,na.rm=T)),color="blue",linetype="solid",size=1) +
     geom_vline(aes(xintercept=mean(LENGTH,na.rm=T)),color="red",linetype="dashed",size=1) +
     geom_hline(aes(yintercept=mean(GC.COUNT/N.COUNT,na.rm=T)),color="red",linetype="dashed",size=1) +
     xlab("Fragment GC content (fraction)") + ylab("Frequency") +
     opts(title=sampleName);
   ggsave(gc.length.file,width=7,height=7,units="in");
  
   save(analysis,file=data.file);
   summary.list = list(L.MEAN=mean(analysis$data$LENGTH),
                               L.MEDIAN=median(analysis$data$LENGTH),
                               L.IQR=analysis$quartiles[4] - analysis$quartiles[2],
                               L.SD=analysis$sd,
                               L.LOW = mean(plot.data$lowEnd),
                               L.HIGH = mean(plot.data$highEnd),
                               NEG.IL.PC=length(which(analysis$data$I.LENGTH < 0))/nrow(analysis$data),
                               GC.MEAN=mean(analysis$data$GC.COUNT/analysis$data$N.COUNT),
                               GC.MEDIAN=median(analysis$data$GC.COUNT/analysis$data$N.COUNT),
                               GC.COR=analysis$qlim.gc.cor.test$estimate,
                               GC.PVALUE=analysis$qlim.gc.cor.test$p.value);
   summary = data.frame(summary.list, row.names=NULL,check.rows=T);
   write.table(summary,file=summary.file,quote=F,row.names=F,sep="\t");
   cat("done\n",file=stderr()); flush(stderr());

} 


referenceSimulation = function(results,simulationCount,excludeMasked=T) {  

  ref.simulation.length = if (simulationCount == nrow(results$data)) results$data$LENGTH else sample(results$data$LENGTH,simulationCount,replace=T);
  ref.simulation.read.length = ceiling(mean(c(results$data$N.FIRST,results$data$N.LAST),na.rm=T));
  ref.simulation.strand = runif(simulationCount) < 0.5;
  ref.simulation.reads.gccount = integer(simulationCount);
  ref.simulation.reads.nucscount = integer(simulationCount);
  ref.simulation.gccount = integer(simulationCount);
  ref.simulation.nucscount = integer(simulationCount);
  ref.simulation.first.gccount = integer(simulationCount);
  ref.simulation.last.gccount = integer(simulationCount);
  ref.simulation.first.nucscount = integer(simulationCount);
  ref.simulation.last.nucscount = integer(simulationCount);
  
  range.widths = sapply(1:length(range.list),simplify=T, function(x) { as.integer(width(range.list[x]))});

  sc.pc = ceiling(simulationCount / 100);
  cat(paste("Performing reference simulation... (",simulationCount," repeats) ",sep=""),file=stderr()); flush(stderr());
  for (i in 1:simulationCount) {
    len = ref.simulation.length[i];
    if (!ref.simulation.strand[i]) len = -len;
    startPos = NULL;
    stopPos = NULL;
    seqName = NULL;
    regionStr = NULL;
    while (is.null(startPos)) {
      candidate = 1 + floor(runif(1) * range.list.total.width);
      range.idx = 1;
      while (range.idx < length(range.list) & candidate > range.widths[range.idx])  {
        candidate = candidate - range.widths[range.idx];
        range.idx = range.idx + 1;
      }
      seqName = names(range.list)[range.idx];
      regionNs = referenceMask.N[[seqName]];
      if (regionNs[candidate] & excludeMasked) next;
      stopPos = candidate + len - 1;
      if (stopPos - ref.simulation.read.length + 1 <= 0 | stopPos > range.widths[range.idx]) next;
      if (regionNs[stopPos] & excludeMasked) next;
      startPos = candidate;
    }
    biasCounts = regionBiasCounts(reference,seqName,startPos,stopPos);
    first.read.counts = regionBiasCounts(reference,seqName,startPos,startPos + ref.simulation.read.length - 1);
    last.read.counts = regionBiasCounts(reference,seqName,stopPos - ref.simulation.read.length + 1,stopPos);
    ref.simulation.gccount[i] = biasCounts[1];
    ref.simulation.nucscount[i] = biasCounts[2];
    ref.simulation.first.gccount[i] = first.read.counts[1];
    ref.simulation.last.gccount[i] = last.read.counts[1];
    ref.simulation.first.nucscount[i] = first.read.counts[2];
    ref.simulation.last.nucscount[i] = last.read.counts[2];
    if (i %% sc.pc == 0) { cat(".",file=stderr()); flush(stderr()) }
  }
  ref.simulation.reads.gccount = ref.simulation.first.gccount + ref.simulation.last.gccount;
  ref.simulation.reads.nucscount = ref.simulation.first.nucscount + ref.simulation.last.nucscount;
  cat("done\n",file=stderr()); flush(stderr());
  res = data.frame(list(GC.COUNT=ref.simulation.gccount,N.COUNT=ref.simulation.nucscount,
                  GC.MAPPED=ref.simulation.reads.gccount,
                  N.MAPPED=ref.simulation.reads.nucscount,
                  GC.FIRST=ref.simulation.first.gccount,
                  N.FIRST=ref.simulation.first.nucscount,
                  GC.LAST=ref.simulation.last.gccount,
                  N.LAST=ref.simulation.last.nucscount));
  res;
}

biasAnalysis = function(results,simulations,excludeMasked=T) {
  ref.simulation.gccount = simulations$GC.COUNT;
  ref.simulation.nucscount = simulations$N.COUNT;
  ref.simulation.reads.gccount = simulations$GC.MAPPED;
  ref.simulation.reads.nucscount = simulations$N.MAPPED;
  ref.simulation.first.gccount = simulations$GC.FIRST;
  ref.simulation.first.nucscount = simulations$N.FIRST;
  ref.simulation.last.gccount = simulations$GC.LAST;
  ref.simulation.last.nucscount = simulations$N.LAST;
  simulationCount = nrow(simulations);
  #FIRST_FACE_AWAY = 1, LAST_FACE_AWAY = 2, FIRST_UNMAPPED = 4, LAST_UNMAPPED = 8,
  #FIRST_LOW_MQ = 16, LAST_LOW_MQ =  32, LARGE_FRAGMENT = 64, CROSS_CHROMOSOME = 128;
  flags = results$data.all$FLAGS;
  flags[!is.na(results$data.all$LENGTH) & results$data.all$LENGTH < results$lowEnd] = 
    bitOr(flags[!is.na(results$data.all$LENGTH) & results$data.all$LENGTH < results$lowEnd],256);
  flags[!is.na(results$data.all$LENGTH) & results$data.all$LENGTH > results$highEnd] = 
    bitOr(flags[!is.na(results$data.all$LENGTH) & results$data.all$LENGTH > results$highEnd],64);
  
  flags.mapped = bitAnd(flags, 4 + 8) == 0;
  flags.both.lowmq = flags.mapped & bitAnd(flags, 16 + 32) == 16 +32;
  flags.one.lowmq = flags.mapped & bitAnd(flags, 16 + 32) != 0 & !flags.both.lowmq; 
  flags.both.unmapped = bitAnd(flags,4) != 0 & bitAnd(flags,8) != 0;
  flags.one.unmapped = !flags.both.unmapped & bitAnd(flags, 4 + 8) != 0;
  flags.sv = flags.mapped & !flags.both.lowmq & !flags.one.lowmq & bitAnd(flags,1 + 2 + 64 + 128 + 256) != 0;
  flags.proper = flags == 0;
  flags.itrans = flags.sv & (flags == 0 | flags == 64) & !is.na(results$data.all$LENGTH) & results$data.all$LENGTH > 10000;
  flags.short = flags.sv & flags == 256;
  flags.long = flags.sv & flags == 64 & !flags.itrans;
  flags.otrans = flags.sv & (bitAnd(flags,128) != 0);
  flags.bfa = flags.sv & bitAnd(flags,1+2) & bitAnd(flags,1 + 2) == 3;
  flags.ofa = flags.sv & (bitAnd(flags,1 + 2) == 1 | bitAnd(flags,1 + 2) == 2);
  
  sequenced.pair.gcratio = ( results$data.all$GC.FIRST + results$data.all$GC.LAST ) / 
                           ( results$data.all$N.FIRST + results$data.all$N.LAST );
  
  mapped.pair.gcratio = sequenced.pair.gcratio[flags.mapped];
  singleton.pair.gcratio = sequenced.pair.gcratio[flags.one.unmapped];
  unmapped.pair.gcratio = sequenced.pair.gcratio[flags.both.unmapped];
  both.lowmq.pair.gcratio = sequenced.pair.gcratio[flags.both.lowmq];
  one.lowmq.pair.gcratio = sequenced.pair.gcratio[flags.one.lowmq];
  proper.pair.gcratio = sequenced.pair.gcratio[flags.proper];
  sv.pair.gcratio = sequenced.pair.gcratio[flags.sv];
  
  reference.fragment.gcratio = ref.simulation.gccount/ ref.simulation.nucscount
  reference.pair.gcratio = ref.simulation.reads.gccount/ref.simulation.reads.nucscount;

  proper.fragment.gcratio = results$data$GC.COUNT[flags.proper]/results$data$N.COUNT[flags.proper];
  
  itrans.pair.gcratio = sequenced.pair.gcratio[flags.itrans];
  otrans.pair.gcratio = sequenced.pair.gcratio[flags.otrans];
  bfa.pair.gcratio = sequenced.pair.gcratio[flags.bfa];
  ofa.pair.gcratio = sequenced.pair.gcratio[flags.ofa];
  short.pair.gcratio = sequenced.pair.gcratio[flags.short];
  long.pair.gcratio = sequenced.pair.gcratio[flags.long];
  
  mapped.read.gcratio = results$data.all$GC.MAPPED[!flags.both.unmapped] / results$data.all$N.MAPPED[!flags.both.unmapped];
  unmapped.read.gcratio = (results$data.all$GC.FIRST[!flags.mapped] + results$data.all$GC.LAST[!flags.mapped] - results$data.all$GC.MAPPED[!flags.mapped]) /
                          (results$data.all$N.FIRST[!flags.mapped]  + results$data.all$N.LAST[!flags.mapped] - results$data.all$N.MAPPED[!flags.mapped]);
  
  gcratio.range = range(c(mapped.read.gcratio,unmapped.read.gcratio));
  gcratio.breaks = seq(floor(gcratio.range[1]),ceiling(gcratio.range[2]),0.01);
  gcratio.hist.breaks = c(gcratio.range[1],(gcratio.breaks[2:(length(gcratio.breaks)-1)] - 0.00001),gcratio.range[2]); 

  both.lowmq.pair.gcratio.hist = hist(both.lowmq.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  one.lowmq.pair.gcratio.hist = hist(one.lowmq.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  mapped.pair.gcratio.hist = hist(mapped.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  singleton.pair.gcratio.hist = hist(singleton.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  unmapped.pair.gcratio.hist = hist(unmapped.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  sv.pair.gcratio.hist = hist(sv.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  proper.pair.gcratio.hist = hist(proper.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  reference.pair.gcratio.hist = hist(reference.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  sequenced.pair.gcratio.hist = hist(sequenced.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  mapped.read.gcratio.hist = hist(mapped.read.gcratio,breaks=gcratio.hist.breaks,plot=F);
  unmapped.read.gcratio.hist = hist(unmapped.read.gcratio,breaks=gcratio.hist.breaks,plot=F);
  
  itrans.pair.gcratio.hist = hist(itrans.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  otrans.pair.gcratio.hist = hist(otrans.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  bfa.pair.gcratio.hist = hist(bfa.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  ofa.pair.gcratio.hist = hist(ofa.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  short.pair.gcratio.hist = hist(short.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  long.pair.gcratio.hist = hist(long.pair.gcratio,breaks=gcratio.hist.breaks,plot=F);
  
  reference.fragment.gcratio.hist = hist(reference.fragment.gcratio,breaks=gcratio.hist.breaks,plot=F);
  proper.fragment.gcratio.hist = hist(proper.fragment.gcratio,breaks=gcratio.hist.breaks,plot=F);
  
  bias.data = list(
      flags=data.frame(
        mapped=flags.mapped,
        both.lowmq=flags.both.lowmq,
        one.lowmq=flags.one.lowmq,
        sv=flags.sv,
        proper=flags.proper,
        itrans=flags.itrans,
        short=flags.short,
        long=flags.long,
        otrans=flags.otrans,
        bfa=flags.bfa,
        ofa=flags.ofa
      ),
      gcratios=gcratio.breaks[1:(length(gcratio.breaks)-1)],
      pair.gcratio.hists = list(
        reference=reference.pair.gcratio.hist,
        sequenced=sequenced.pair.gcratio.hist,
        mapped=mapped.pair.gcratio.hist,
        unmapped=unmapped.pair.gcratio.hist,
        proper=proper.pair.gcratio.hist,
        singleton=singleton.pair.gcratio.hist,
        sv=sv.pair.gcratio.hist,
        itrans=itrans.pair.gcratio.hist,
        otrans=otrans.pair.gcratio.hist,
        bfa=bfa.pair.gcratio.hist,
        ofa=ofa.pair.gcratio.hist,
        short=short.pair.gcratio.hist,
        blowmq=both.lowmq.pair.gcratio.hist,
        olowmq=one.lowmq.pair.gcratio.hist,
        long=long.pair.gcratio.hist),
      fragment.gcratio.hists = list(
        reference=reference.fragment.gcratio.hist,
        proper=proper.fragment.gcratio.hist),
      pair.gcratios=list(
        reference=reference.pair.gcratio,
        sequenced=sequenced.pair.gcratio,
        mapped=mapped.pair.gcratio,
        unmapped=unmapped.pair.gcratio,
        blowmq=both.lowmq.pair.gcratio,
        olowmq=one.lowmq.pair.gcratio,
        proper=proper.pair.gcratio,
        singleton=singleton.pair.gcratio,
        sv=sv.pair.gcratio,
        itrans=itrans.pair.gcratio,
        otrans=otrans.pair.gcratio,
        bfa=bfa.pair.gcratio,
        ofa=ofa.pair.gcratio,
        short=short.pair.gcratio,
        long=long.pair.gcratio),
      fragment.gcratios=list(
        reference=reference.fragment.gcratio,
        proper=proper.fragment.gcratio)
  );
}

writeBiasResults = function(bias,outDir,sampleName="") {
  # Absolute bias, Reference vs Sample sequencing and mapped based on fragment composition.
  bias.plot.file = paste(outDir,"bias.svg",sep="/");
  fragments.density.plot.file = paste(outDir,"fragments.density.svg",sep="/");
  bias.ratios.log.plot.file = paste(outDir,"bias.ratios.log.svg",sep="/");
  bias.ratios.plot.file = paste(outDir,"bias.ratios.svg",sep="/");
  bias.data.file = paste(outDir,"bias.RData",sep="/");
  bias.summary.file = paste(outDir,"bias.tsv",sep="/");
  bias.seq.vs.ref.plot.file = paste(outDir,"reference-vs-sequenced.fractions.svg",sep="/");
  bias.seq.vs.ref.f.plot.file = paste(outDir,"reference-vs-sequenced-fragments.fractions.svg",sep="/");
  bias.seq.vs.proper.plot.file = paste(outDir,"sequenced.fractions.svg",sep="/");
  bias.sv.types.plot.file = paste(outDir,"sv.fractions.svg",sep="/");
  
  proper.fragment.density = bias$fragment.gcratio.hists$proper$density;
  reference.fragment.density = bias$fragment.gcratio.hists$reference$density;
  
  cat("Plotting simulations bias...",file=stderr()); flush(stderr());
  
  fragments.data = data.frame(gc=bias$gcratios,reference=reference.fragment.density,proper=proper.fragment.density);
  fragments.data$ratio = fragments.data$proper/fragments.data$reference;
  fragments.data$ratio[!is.finite(fragments.data$ratio)] = NA;
  fragments.data$reference = (fragments.data$reference / max(fragments.data$reference)) * 3;
  fragments.data$proper = (fragments.data$proper / max(fragments.data$proper)) * 3;
  
  ggplot(data=fragments.data,aes(x=gc)) + xlab("GC fragment composition (%)") + ylab("Density ratio (Sample / Reference)")  +
    stat_smooth(aes(y=ratio),span=0.3,na.rm=T,method="loess",se=T) + geom_point(aes(y=ratio),na.rm=T) +
    stat_smooth(aes(y=reference),color="red",span=0.1,method="loess",na.rm=T,se=F) +
    stat_smooth(aes(y=proper),color="green",span=0.1,method="loess",na.rm=T,se=F) +
    xlab("Fragment CG Content (fraction)") + ylab("Relative density; Ratio proper fragment vs reference") +
    opts(title=sampleName);
  ggsave(bias.plot.file,width=7,height=7,units="in");
  save(bias,file=bias.data.file);

  densities = list(
    reference=density(bias$pair.gcratios$reference,na.rm=T),
    sequenced=density(bias$pair.gcratios$sequenced,na.rm=T),
    mapped=density(bias$pair.gcratios$mapped,na.rm=T),
    proper=density(bias$pair.gcratios$proper,na.rm=T),
    unmapped=density(bias$pair.gcratios$unmapped,na.rm=T),
    sv=density(bias$pair.gcratios$sv,na.rm=T));
    
  ggReadsGCDensityRatioPlot(densities,control="reference",log=T) +
    xlab("Read pair GC content (fraction)") + ylab("Density ratio versus the reference (log)") +
    opts(titel=sampleName) +  
    ggsave(bias.ratios.log.plot.file,width=10,height=7,units="in");

  ggReadsGCDensityRatioPlot(densities,control="reference",ylim=c(0,4)) +
    xlab("Read pair GC content (fraction)") + ylab("Density ratio versus the reference") +
    opts(titel=sampleName) +  
  ggsave(bias.ratios.plot.file,width=10,height=7,units="in");

  ref.vs.seq.ffractions = data.frame(
    proper=bias$fragment.gcratio.hists$proper$counts,
    reference=bias$fragment.gcratio.hists$reference$counts);
  
  ref.vs.seq.ffractions.colSum = colSums(ref.vs.seq.ffractions);
  ref.vs.seq.ffractions.divisor = ref.vs.seq.ffractions.colSum / min(ref.vs.seq.ffractions.colSum);
  ref.vs.seq.ffractions = ref.vs.seq.ffractions / matrix(rep(ref.vs.seq.ffractions.divisor,nrow(ref.vs.seq.ffractions)),nrow=nrow(ref.vs.seq.ffractions),byrow=T);
  ggReadsGCFrequencyPartitionPlot(ref.vs.seq.ffractions,bias$gcratios,minSum=100,smooth="proper") +
    xlab("Fragment GC content (fraction)") + ylab("Count fraction Proper vs reference") + 
    opts(title=sampleName); 
  ggsave(bias.seq.vs.ref.f.plot.file,width=10,height=7,units="in");
  
  ref.vs.seq.fractions = data.frame(
          sequenced=bias$pair.gcratio.hists$sequenced$counts,
          reference=bias$pair.gcratio.hists$reference$counts);
  
  ref.vs.seq.fractions.colSum = colSums(ref.vs.seq.fractions);
  ref.vs.seq.fractions.divisor = ref.vs.seq.fractions.colSum / min(ref.vs.seq.fractions.colSum);
  ref.vs.seq.fractions = ref.vs.seq.fractions / matrix(rep(ref.vs.seq.fractions.divisor,nrow(ref.vs.seq.fractions)),nrow=nrow(ref.vs.seq.fractions),byrow=T);
  ggReadsGCFrequencyPartitionPlot(ref.vs.seq.fractions,bias$gcratios,minSum=100,smooth="sequenced")  +
    xlab("Read pair GC content (fraction)") + ylab("Count fraction Proper vs reference") + 
    opts(title=sampleName); 
  ggsave(bias.seq.vs.ref.plot.file,width=10,height=7,units="in");

  seqs.fractions = data.frame(
    proper=bias$pair.gcratio.hists$proper$counts,
    one.lowmq=bias$pair.gcratio.hist$olowmq$counts,
    both.lowmq=bias$pair.gcratio.hists$blowmq$counts,
    unmapped=bias$pair.gcratio.hists$unmapped$counts,
    singleton=bias$pair.gcratio.hists$singleton$counts,
    sv=bias$pair.gcratio.hists$sv$counts);
  ggReadsGCFrequencyPartitionPlot(seqs.fractions,bias$gcratios,minSum=100)  +
    xlab("Mapped read-pair GC content (fraction)") + ylab("Count fraction between mapped pair categories") + 
    opts(title=sampleName);
  ggsave(bias.seq.vs.proper.plot.file,width=10,height=7,units="in");
  
  sv.fractions = data.frame(
    long=bias$pair.gcratio.hists$long$counts,
    short=bias$pair.gcratio.hists$short$counts,
    both.face.away=bias$pair.gcratio.hists$bfa$counts,
    one.face.away=bias$pair.gcratio.hists$ofa$counts,
    same.chr.trans=bias$pair.gcratio.hists$itrans$counts,
    cross.chr.trans=bias$pair.gcratio.hists$otrans$counts
  );
  
  ggReadsGCFrequencyPartitionPlot(sv.fractions,bias$gcratios,minSum=100)  +
    xlab("Mapped read-pair GC content (fraction)") + ylab("Count fraction between mapped pair categories") + 
    opts(title=sampleName);
  ggsave(bias.sv.types.plot.file,width=10,height=7,units="in");
  
  bias.data = data.frame(gc=bias$gcratios,
                         fragments.reference=fragments.data$reference,
                         fragments.proper=fragments.data$proper,
                         fragments.proper.vs.reference.ratio=fragments.data$ratio);
                              
  write.table(bias.data,file=bias.summary.file,quote=F,row.names=F,sep="\t");
  cat("done\n",file=stderr()); flush(stderr());
}

gcCorrection = function(results, bin.size=1000,outDir,sampleName='') {
  if (!require("GCcorrect")) stop("GCcorect is not present");
  basic.plots.file = paste(outDir,"GCcorrection.basic.pdf",sep="/");
  comp.plots.file = paste(outDir,"GCcorrection.comp.svg",sep="/");
  gcwindow.file.prefix = paste(outDir,"GCcorrection.gcwindow",sep="/");
  tvs.plot.file = paste(outDir,"GCcorrection.tvs.svg",sep="/");
  cmeans.plot.file = paste(outDir,"GCcorrection.cmeans.svg",sep="/");
  norm.plot.file = paste(outDir,"binnorm.svg",sep="/");
  norm.hist.plot.file = paste(outDir,"binhist.svg",sep="/");
  norm.gc.plot.file = paste(outDir,"bingc.svg",sep="/");
  correction.data.file = paste(outDir,"correction.RData",sep="/");
  isrep = c(referenceMask.N , recursive=T);
  mask.chr = c(sapply(reference,simplify="array",function(x) { chars = strsplit(as.character(x),""); }),recursive=T);
  chrline = as.integer(factor(mask.chr,levels=c("A","T","C","G"))) - 1;
  chr = prepareChrom("/dev/null",clean=T,repeats=isrep,reference=chrline,savefile=F);
  so.far = 0;
  concat.offsets = integer(sum(width(reference)));
  concat.offsets[1:length(concat.offsets)] = NA;
  so.far.offset = 0;
  chroffsets = integer(length(reference));
  for (i in 1:length(reference)) {
    full.range = c(1,width(reference[i]));
    chr.name = names(reference[i]);
    full.irange = IRanges(start=1,end=width(reference[i]));
    full.width = sum(sum(width(full.irange)));
    interval.irange = range.list[[chr.name]];
    diffs = if (is.null(interval.irange)) full.irange else setdiff(full.irange,interval.irange);
    if (length(diffs) > 0) {
      for (j in 1:length(diffs)) {
        diff.range = c(start(diffs[j]),end(diffs[j])) + so.far;
        isrep[diff.range[1]:diff.range[2]] = T;
        diff.width = diff.range[2] - diff.range[1] + 1;
        concat.offsets[diff.range[1]:diff.range[2]] = c((so.far.offset + 1):(so.far.offset + diff.width));
        so.far.offset = so.far.offset + diff.width;
      }
    }
    chroffsets[i] = so.far;
    so.far = so.far + full.width;
  }
  names(chroffsets) = names(reference);
  dat = results$data;
  for (i in 1:length(chroffsets)) {
    chr.rows = which(dat$CHROM == names(chroffsets)[i]);
    dat$START[chr.rows] = dat$START[chr.rows] + chroffsets[i]; 
    dat$END[chr.rows] = dat$END[chr.rows] + chroffsets[i]; 
  }
  
  # custom loess #
  bin.count = ceiling(length(isrep)/bin.size);
  dat.is.for.strand = dat$STRAND == "+";
  bin.num = integer(nrow(dat));
  bin.num[dat.is.for.strand] = floor((dat$START[dat.is.for.strand]-1)/bin.size) + 1;
  bin.num[!dat.is.for.strand] = floor((dat$END[!dat.is.for.strand]-1)/bin.size) + 1;
  bin.n.total = sapply(1:bin.count,function(i) { length(which(!isrep[(1+(i-1)*bin.size):(i*bin.size)])) });
  bin.starts = integer(bin.count);
  bin.gc.count = integer(bin.count);
  bin.n.count = integer(bin.count);
  bin.gc.bias.sum = numeric(bin.count);
  bin.gc.bias.count = integer(bin.count);
  for (i in 1:nrow(dat)) {
    bn = bin.num[i];
    bin.starts[bn] = bin.starts[bn] + 1;
    bin.gc.count[bn]  = bin.gc.count[bn] + dat$GC.COUNT[i];
    bin.n.count[bn] = bin.n.count[bn] + dat$N.COUNT[i];
    if (dat$N.COUNT[i] > 0) {
      bin.gc.bias.sum[bn] = bin.gc.bias.sum[bn] + (dat$GC.COUNT[i] / dat$N.COUNT[i])
      bin.gc.bias.count[bn] = bin.gc.bias.count[bn] + 1;
    }
  }
  
  bins = data.frame(NUM=1:bin.count,SIZE=bin.n.total,STARTS=bin.starts,GC.COUNT=bin.gc.count,N.COUNT=bin.n.count,
                    GC.BIAS=bin.gc.bias.sum / bin.gc.bias.count);
  bins$EXCLUDE = is.na(bins$GC.BIAS) | bins$SIZE < 0.5 * bin.size;
  bins.loess = loess((STARTS / SIZE) ~ GC.BIAS,data=bins,subset=!EXCLUDE,span=0.3);
  bins$PRED = predict(bins.loess,bins$GC.BIAS)
  bins$CORRECTION = bins$STARTS / (bins$PRED * bins$SIZE);  

  #Gccorrect
  dat1_for = as.matrix(dat[dat$STRAND == "+",2:3]);
  dat1_rev = as.matrix(dat[dat$STRAND == "-",2:3]);
  
  
  dat1 = prepareReads (dat1_for, dat1_rev,chr_len = length(chr$chrline), 0,76,'Plasodium_falciparum');  
  pdf(file=basic.plots.file,width=9,height=12);
    basicPlots(dat1,types= c(1,2,3,4), binsize = "10K");
  dev.off();
  #  compDataRefPlots(chr,dat1,types = c(1,2,3,4));
  #dev.off();

  useSamp = !logical(length(chr$isrep));
  margin = 5;
  sampCh1 = sampleChrom(chr,dat1,n=simulationCount,margin = margin,len_range = c(), memoryopt = TRUE,useSamp=useSamp); 
  begdata = makeGCLens(chr$isgc,dat1$forw,sampline = sampCh1$singleLocSamp,minlen = 8,maxlens=500,margin=margin,max_frag_for_loc=10)
  tvs = scoreGCLens(begdata, maxlen=500, minlen =  8,scale= T);
  svg(file=tvs.plot.file,width=5,height=4);
    plotGCLens(tvs,lw=2,lt=2)
  dev.off();
  best_size = which.max(tvs)-1;
  #best_size = results$median;
  gcsize = best_size;
  gcline = prepareGCChrom(chr,gcsize,filename = gcwindow.file.prefix);
  cMeans= getCondMean(gcline[sampCh1$singleLocSamp+margin],sampCh1$forsamped,cutoff = 4,jump = 6,norm = FALSE);
  cMeans_rev= getCondMean(gcline[pmax(1,sampCh1$singleLocSamp+dat1$readlen-1-margin-gcsize)],sampCh1$revsamped,cutoff = 4,jump = 6,norm = FALSE);
  svg(cmeans.plot.file,width=7,height=5);
    par(mfcol = c(1,2))
    plotCondMean(cMean = cMeans,ci = TRUE,normRange = gcsize,meanLine=TRUE,lt = 1,col=4)
    plotCondMean(cMean = cMeans_rev,ci = TRUE,normRange = gcsize,meanLine=TRUE,lt = 1,col=4)
  dev.off();
  cMeans_tog = sumCondMeans(cMeans ,cMeans_rev);
  forpreds = predictLine(cMeans_tog,gcline,gcsize,margin,chr$isrep,strand = "F",readlen=dat1$readlen);
  revpreds = predictLine(cMeans_tog,gcline,gcsize,margin,chr$isrep,strand = "R",readlen=dat1$readlen);
  
  pred1K = meanLine(forpreds$preds,bin.size)*bin.size;
  pred1K_rev = meanLine(revpreds$preds,bin.size)*bin.size;
  
  pred1K = meanLine(forpreds$preds,bin.size)*bin.size;
  pred1K_rev= meanLine(revpreds$preds,bin.size)*bin.size;
  
  chr$gc1K = meanLine(chr$isgc,bin.size)
  chr$map1K = meanLine(chr$isrep,bin.size)
  dat1$for1K = binReads(dat1$forw[,1],bin.size,dat1$chr_len)
  dat1$rev1K = binReads(dat1$reve[,1],bin.size,dat1$chr_len)
  
  norm1K = (dat1$for1K+dat1$rev1K)/(pred1K+pred1K_rev+0.01);
  unnorm1K = (dat1$for1K+dat1$rev1K)/median(dat1$for1K+dat1$rev1K, na.rm=TRUE);

  dat1$both1K = dat1$for1K + dat1$rev1K;
  loessModel = loess(dat1$both1K[dat1$both1K > 0] ~ chr$gc1K[dat1$both1K > 0],span=0.2);
  loessPred = predict(loessModel,chr$gc1K);
  
  dat1$corrected = dat1$both1K / loessPred;

  svg(file=norm.plot.file,width=36,height=16);
    layout(c(1:4));
     normPlot(unnorm1K,chr$gc1K,main=paste(sampleName," Unormalized ",bin.size,"bp bins",sep=""));
     normPlot(norm1K,chr$gc1K,main=paste(sampleName," GCcorrect normalized ",bin.size,"bp bins",sep=""));
     normPlot(dat1$corrected,chr$gc1K,main=paste(sampleName," Loess by bin GC normalized ",bin.size,"bp bins",sep=""));
     normPlot(bins$CORRECTION,bins$GC.BIAS,main=paste(sampleName," Loess by fragment GC average ",bin.size,"bp bins",sep=""));
  dev.off();
  
  norm.hist.plot.file = paste(outDir,"binhist.svg",sep="/");
  norm.gc.plot.file = paste(outDir,"bingc.svg",sep="/");
  correction.data.file = paste(outDir,"correction.RData",sep="/");
  
  svg(file=norm.gc.plot.file,width=10,height=15);
    ylim=c(0,3);
    xlab="GC content";
    ylab="scaled fragment starts (mean to 1)";
    layout(matrix(c(1:6),byrow=T,nrow=3,ncol=2));
    plot(chr$gc1K,unnorm1K,xlab=xlab,ylab=ylab,ylim=ylim,pch='.',main=paste(sampleName," Unormalized ",bin.size,"bp bin fragment starts vs GC content",sep=""));
    points(chr$gc1K,(pred1K+pred1K_rev)/median(pred1K+pred1K_rev),col="blue",pch=".");
    plot(chr$gc1K,norm1K,xlab=xlab,ylab=ylab,ylim=ylim,pch='.',main=paste(sampleName,"GCcorrect normalized ",bin.size,"bp bin fragment starts vs GC content",sep=""));
    plot(chr$gc1K,unnorm1K,xlab=xlab,ylab=ylab,ylim=ylim,pch='.',main=paste(sampleName," Unormalized ",bin.size,"bp bin fragment starts vs GC content",sep=""));
    points(chr$gc1K,loessPred/median(loessPred),col="red",pch=".");
    plot(chr$gc1K,dat1$corrected,xlab=xlab,ylab=ylab,ylim=ylim,pch='.',main=paste(sampleName," Loess by bin GC normalized ",bin.size,"bp bins",sep=""));
    
    plot(bins$GC.BIAS,(bins$STARTS / bins$SIZE)/median(bins$STARTS/bins$SIZE,na.rm=T),xlab=xlab,ylab=ylab,ylim=ylim,pch='.',main=paste(sampleName," Uncorrected number of starts vs fragment GC average ",bin.size,"bp bin fragment starts vs GC content",sep=""));
    points(bins$GC.BIAS,bins$PRED/median(bins$PRED,na.rm=T),col="orange",pch=".");
    plot(bins$GC.BIAS,bins$CORRECTION,xlab=xlab,ylab=ylab,ylim=ylim,pch='.',main=paste(sampleName," Loess by fragment GC average ",bin.size,"bp bins",sep=""));
  dev.off();
}

normPlot = function(y,gc,main="",pch='.') {
plot(y,cex = 0.5, pch='.',main=main,ylim=c(0,3),col="white");
segments(-length(y),1,0,1,lwd=1);
segments(length(y),1,length(y)*2,lwd=1)
abline(h=2);
x = y; x[gc > 0.1] = NA; points(x,cex = 0.5, pch=pch,ylim=c(0,3),col="orange");
x = y; x[gc < 0.1 | gc > 0.2] = NA; points(x,cex = 0.5, pch=pch,ylim=c(0,3),col="red");
x = y; x[gc < 0.2 | gc > 0.3] = NA; points(x,cex = 0.5, pch=pch,ylim=c(0,3),col="green");
x = y; x[gc < 0.3] = NA; points(x,cex = 0.5, pch=pch,ylim=c(0,3),col="blue");
}


ggReadsGCDensityPlot = function(densities,names=NULL,colours=rainbow(length(densities)+1),size=2) {
  if (is.null(names)) names = names(densities);
  if (is.null(names)) names = c(1:length(densities));
  colours = rep(colours,length(densities));
  densities.length = length(densities);
  dframes = lapply(c(1:densities.length),function(i) {
    d = densities[[i]];
    l = list(x=d$x,y=d$y,name=names[i]);
    data.frame(l);
  });
  
  all.dframe = do.call(rbind,dframes);
  plot = ggplot(data=all.dframe,aes(x=x,y=y,group=name)) + geom_line(aes(color=name),size=2) + 
    scale_colour_manual("",breaks=names,values=colours) + opts(legend.text = theme_text(size=20));
  plot
}

ggReadsGCDensityRatioPlot = function(densities,control,names=NULL,points=100,log=F,ylim=NULL,...) {
  if (is.null(names)) names = names(densities);
  if (is.null(names)) names = c(1:length(densities));
  if (is.null(control)) stop("you must specify a control density index or name");
  dframes.0 = lapply(densities,function(d) {
    res = data.frame(list(x=d$x,y=d$y));
  });
  control.name = if (is.integer(control)) names[control] else control;
  control.idx = if (is.integer(control)) control else match(control,names);
  control.df = dframes.0[[control.idx]];
  dframes = lapply(dframes.0, function(d) { 
    print(range(d[,1]));
    print(range(control.df[,1]))
    res = lineRatio::lineRatio.density(d,control.df);
    colnames(res) = c("x","y"); 
    as.data.frame(res) 
  });
  m = ggReadsGCDensityPlot(densities=dframes,names=names) + ylab("Ratio versus");
  if (!is.null(ylim)) m = m + coord_cartesian(ylim=ylim);
  if (log) m = m + scale_y_log10(breaks=c(0,0.1,0.25,0.5,1,2,4,10));
  m
}

ggReadsGCFrequencyPartitionPlot = function(freqs,breaks,colours=rainbow(nrow(shares)+1),minSum="auto",smooth=NULL,loess=NULL) {
  if (is.null(freqs)) stop("you must specify a frequency list");
  if (is.null(breaks)) stop("you must sepcify the breaks");
  len = length(breaks);
  if (nrow(freqs) != len) stop("number of rows in freqs does not match the number of breaks");
  total = rowSums(freqs,na.rm=T);
  bar.width = integer(length(total)) + 1;
  if (!is.null(minSum)) {
    if (minSum == "auto") minSum = max(100,sqrt(sum(total)));
    accu = 0;
    curr.group = 1;
    group = integer(length(breaks));
    for (i in 1:length(group)) {
      group[i] = curr.group;
      accu  = accu + total[i];
      if (i != 1) { if (total[i] > total[i-1] * 10) next; if (total[i] < total[i-1] * 0.1) next };
      if (accu >= minSum) { accu = 0; curr.group = curr.group + 1; }
    }
    if (accu < (minSum * 0.5)) { curr.group = curr.group - 1; group[group == curr.group +1] = curr.group; }
    bar.width = sapply(c(1:max(group)),simplify=T,function(g) { 
      low = if (g == 1) 0 else (max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) /2
      high = if (g == max(group)) max(breaks) + (-max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) / 2 else (max(breaks[which(group==g)]) + min(breaks[which(group==g+1)])) /2
      high - low
    });
    breaks = sapply(c(1:max(group)),simplify=T,function(g) { 
      low = if (g == 1) 0 else (max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) /2
      high = if (g == max(group)) max(breaks) + (-max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) / 2 else (max(breaks[which(group==g)]) + min(breaks[which(group==g+1)])) /2
      (high + low) /2
    }); 
    if (curr.group < len) {
      freqs = data.frame(sapply(names(freqs),USE.NAMES=T,function(x) {
        list(sapply(c(1:max(group)),simplify=T,function(g) { sum(freqs[group == g,x]) }));
      }));
      total = rowSums(freqs,na.rm=T);
    }
    
  }
  dframes = lapply(c(1:ncol(freqs)),function(i) {
    f = freqs[,i] / total;
    l = list(gc=breaks,bw=bar.width,frac=f,name=names(freqs)[i]);
    data.frame(l);
  });
  data = do.call(rbind,dframes);
  print(data);
  res = ggplot(data=data,aes(x=gc,y=frac,fill=name)) + ylab("Fraction") + xlab("GC % in read-pair") + 
    geom_bar(aes(width=bw),stat="identity") #geom_bar(aes(width=bw),stat="identity")
  if (!is.null(smooth)) res =
    res + stat_smooth(data=data[data$name == smooth,],aes(x=gc,y=frac),method="loess",span=0.7) + geom_point(data=data[data$name == smooth,],aes(x=gc,y=frac))
  #      res + stat_smooth(data=data[data$name == smooth,],aes(x=gc,y=frac),method="loess",span=0.3) + geom_points(data=data[data$name == smooth,],aes(x=gc,y=frac))
  res
}


#require("bitops",quietly=T);
#require("ggplot2",quietly=T);
#require("mclust",quitetly=T);
#require("rtracklayer",quietly=T);
#require("IRanges",quietly=T);
#require("Biostrings",quietly=T);

