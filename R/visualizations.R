
ggReadsGCDensityPlot = function(densities,names=NULL,colours=rainbow(length(densities)+1),size=2) {
  if (is.null(names)) names = names(densities);
  if (is.null(names)) names = c(1:length(densities));
  colours = rep(colours,length(densities));
  densities.length = length(densities);
  dframes = lapply(c(1:densities.length),function(i) {
    d = densities[[i]];
    l = list('x'=d$x,'y'=d$y,'name'=names[i]);
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

ggReadsGCFrequencyPartitionPlot = function(freqs,breaks,colours=rainbow(ncol(freqs)+1),
                                           minSum="auto",smooth=NULL,loess=NULL) {
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
      low = if (g == 1) 0 else (max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) / 2
      high = if (g == max(group)) { 
                max(breaks) + (-max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) / 2  
             } else { 
               (max(breaks[which(group==g)]) + min(breaks[which(group==g+1)])) /2 
      };
      high - low
    });
    breaks = sapply(c(1:max(group)),simplify=T,function(g) { 
      low = if (g == 1) 0 else (max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) /2
      high = if (g == max(group)) { max(breaks) + (-max(breaks[which(group==g-1)]) + min(breaks[which(group==g)])) / 2 
           } else { (max(breaks[which(group==g)]) + min(breaks[which(group==g+1)])) /2 }
      (high + low) / 2
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
    l = list('gc'=breaks,'bw'=bar.width,'frac'=f,'name'=names(freqs)[i]);
    data.frame(l);
  });
  data = do.call(rbind,dframes);
  print(data);
  res = ggplot(data=data,aes(x=gc,y=frac,fill=name)) + ylab("Fraction") + xlab("GC % in read-pair") + 
    geom_bar(aes(width=bw),stat="identity");
  if (!is.null(smooth)) res =
    res + stat_smooth(data=data[data$name == smooth,],aes(x=gc,y=frac),method="loess",span=0.7) + geom_point(data=data[data$name == smooth,],aes(x=gc,y=frac))
  res
}
