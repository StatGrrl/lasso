# Plot sample curves and average
plot.path <- function(data, xname, ynames, col=rep("black",length(ynames)), 
                      show.min=rep(T,length(ynames)), 
                      alpha=255, tint=-0.75, ave.only=F, legend=T, 
                      topx=NULL, toplab=NULL, logx=F, ylim=NULL, 
                      lty=rep(1,length(ynames)), lab=c("",""), 
                      title=NULL, ...) {
  # INPUT:
  # data........list containing plot data
  #             each list item is a data frame or matrix of one sample
  #             x and y are named columns in the data frame or matrix
  #             x vector must be identical in each sample
  # xnames......x variable name
  # ynames......vector of y variable names
  # col.........vector of color names for y variables
  # show.min....vector of logicals to plot point at minimum of y variables
  # alpha.......transparency of sample curves (default=none)
  #             integer between 0 (complete) and 255 (none)
  # tint........amount to tint sample curves (default=75% lighter)
  #             float between -1 (white) and 1 (black)
  # ave.only....logical, print only average curves? 
  # legend......logical, display legend? default TRUE
  # topx........to plot secondary x axis on top, vector of same length as x
  # toplab......label for secondary x axis
  # logx........logical, plot the negative log of x variable, default FALSE
  # ylim........vector specified minimum and maximum for y variables
  #             by default, it is calculated from the data
  # lty.........vector of line types for average curves
  # lab.........vector (xlab, ylab) (optional)
  # title.......plot title (optional)
  # ... other arguments to function plot
  #
  # OUTPUT:
  # Plot of sample curves and their averages
  
  source("color.R")

  # graphics parameters
  par(mar=c(3,3,1,1)+0.1, mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  if (legend) {
    legnames <- gsub("."," ", ynames, fixed=T)
    legmax <- max(strwidth(paste0("cccc",legnames), "inches"))
    par(mai=par("mai")+c(0,0,0,legmax))
  }  
  if (as.character(lab[1])=="") par(mar=par("mar")-c(1.5,0,0,0))
  if (as.character(lab[2])=="") par(mar=par("mar")-c(0,1.5,0,0))
  if (!is.null(topx)) par(mar=par("mar")+c(0,0,0.5,0))
  if (!is.null(title)) par(mar=par("mar")+c(0,0,1.5,0))
  
  # get averages
  data2 <- do.call("rbind", data)
  x <- data[[1]][, xname]
  if (logx) x <- -log(x)
  n <- length(x)
  ynum <- length(ynames)
  ave <- data.frame(x, matrix(0, n, ynum))
  names(ave) <- c(xname, ynames)
  for (i in 1:ynum) {
    ave[,i+1] <- tapply(data2[,ynames[i]], data2[,xname], mean)
  }
  
  # y limits
  if (!is.null(ylim)) {
    ymin <- ylim[1]
    ymax <- ylim[2]
  } else if (ave.only) {
      ymin <- floor(min(ave[,ynames]))
      ymax <- ceiling(max(ave[,ynames]))
  } else {
    ymin <- floor(min(data2[,ynames]))
    ymax <- ceiling(max(data2[,ynames]))  
  }
  xlab <- if (is.null(topx)) lab[1] else ""
  plot(x=x, y=rep(1,n), type="n", 
       xlab=xlab, ylab=lab[2], ylim=c(ymin,ymax), ...)
  if (!is.null(topx)) 
    mtext(lab[1], side=1, line=0.25, at=par("usr")[2], las=1)
  
  # sample curves
  if (!ave.only) {
    N <- length(data)
    col.tint <- rep(0, ynum)
    for (i in 1:ynum) {
      col.tint[i] <- color.tint(col[i], alpha, tint)
      for (j in 1:N) {
        lines(x, data[[j]][,ynames[i]], 
              col=col.tint[i]) # sample curves lighter
      }
    }
  }
  # average curves
  for (i in 1:ynum) {
    lines(x, ave[,ynames[i]], col=col[i], lwd=2, lty=lty[i])
    if (show.min[i]) points(x[which.min(ave[,ynames[i]])], 
                            min(ave[,ynames[i]]), col=col[i], pch=19)
  }
  # legend
  if (legend) {
    xpos <- par("usr")[2]
    ypos <- par("usr")[4]
    legend(xpos, ypos, legend=legnames, cex=0.8, 
           col=col, lwd=2, lty=lty, bty="n", xpd=T)
  }
  
  # secondary x axis
  if (!is.null(topx)) {
    x2 <- topx
    x1 <- x
    if (-Inf %in% x) x1 <- x[-which(x==-Inf)]
    if (Inf %in% x) x1 <- x[-which(x==Inf)]
    tick <- axisTicks(c(min(x1),max(x1)), log=F)
    x1 <- sort(x1)
    ind <- match(x1[findInterval(tick,x1)], x)
    axis(side=3, at=x[ind], labels=round(x2[ind]))
    mtext(toplab, side=3, line=0.25, adj=1 )
  }

  # title
  if (!is.null(title)) {
    line <- if (!is.null(topx)) 1.5 else 1
    title(title, line=line)
  }}


# Plot lines
plot.lines <- function(data, xname, ynames, col=rep("black",length(ynames)), 
                       show.min=rep(T,length(ynames)), 
                       lwd=rep(2,length(ynames)), lty=rep(1,length(ynames)), 
                       pch=NULL, type="l", hline=NULL, hlty=NULL, 
                       vline=NULL, vlty=NULL, vlab=NULL, topx=NULL,
                       legend=T, logx=F, ylim=NULL, lab=c("",""), 
                       title=NULL, scale=NULL, ...) {
  # INPUT:
  # data........data frame or matrix where
  #             x and y are named columns
  # xnames......x variable name
  # ynames......vector of y variable names
  # col.........vector of color names for y variables
  # show.min....vector of logicals to plot point at minimum of y variables
  # type........l - lines (default) / b - lines and points
  # lwd.........vector of line widths for y variables, default 2 for all
  # lty.........vector of line types for y variables, default 1 for all
  # pch.........vector of point types for y variables (only if type="b")
  # hline.......vector of points at which to draw horizontal lines
  # hlty........vector of line types for hline
  # vline.......vector of points at which to draw vertical lines
  # vlty........vector of line types for vline
  # vlab........vector of labels for vline
  # topx........to plot secondary x axis on top, provide c(name, label)
  # ylim........specified minimum and maximum for y variables
  # legend......display legend? default TRUE
  # logx........plot the negative log of x variable, default FALSE
  # lab.........vector (xlab, ylab), default no labels
  # title.......plot title, default no title
  # scale.......vector to scale y variables
  # ... other arguments to function plot
  #
  # OUTPUT:
  # Plot lines of y variables

  # graphics parameters
  par(mar=c(3,3,1,1)+0.1, mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  if (legend) {
    legnames <- gsub("."," ", ynames, fixed=T)
    legmax <- max(strwidth(paste0("cccc",legnames), "inches"))
    par(mai=par("mai")+c(0,0,0,legmax))
  }  
  if (as.character(lab[1])=="") par(mar=par("mar")-c(1.5,0,0,0))
  if (as.character(lab[2])=="") par(mar=par("mar")-c(0,1.5,0,0))
  if (!is.null(topx)) par(mar=par("mar")+c(0,0,0.5,0))
  if (!is.null(title)) par(mar=par("mar")+c(0,0,1.5,0))
    
  yvars <- data[,ynames]
  if (length(ynames) == 1) yvars <- matrix(yvars,ncol=1)
  if (!is.null(scale)) scale(yvars, center=F, scale=scale)
  
  ynum <- length(ynames)
  if (!is.null(ylim)) {
    ymin <- ylim[1]
    ymax <- ylim[2]
  } else {
    ymin <- min(yvars)
    ymax <- max(yvars)
  }
  # plot lines
  x <- data[,xname]
  if (logx) x <- -log(x)
  xmax <- max(x)
  n <- length(x)
  xlab <- if (is.null(topx)) lab[1] else ""
  plot(x=x, y=rep(1,n), type="n", 
       xlab=xlab, ylab=lab[2], ylim=c(ymin,ymax), ...)
  if (!is.null(topx)) 
    mtext(lab[1], side=1, line=0.25, at=par("usr")[2], las=1)
  for (i in 1:ynum) {
    lines(x, yvars[,i], type=type, col=col[i], 
          lwd=lwd[i], lty=lty[i], pch=pch[i])
    if (show.min[i]) points(x[which.min(yvars[,i])], 
                            min(yvars[,i]), col=col[i], pch=19)
  }
  # horizontal and vertical lines
  if (!is.null(hline)) 
    for (i in 1:length(hline)) abline(h=hline[i],lty=hlty[i])
  if (!is.null(vline)) 
    for (i in 1:length(vline)) {
      abline(v=vline[i],lty=vlty[i])
      if (!is.null(vlab))
        text(vline[i]*1.01, ymax, labels=vlab[i], adj=0, cex=0.8)
    }
  
  # legend
  if (legend) {
    xpos <- par("usr")[2]
    ypos <- par("usr")[4]
    legend(xpos, ypos, legend=legnames, cex=0.8, pch=pch,
           col=col, lwd=lwd, lty=lty, bty="n", xpd=T)    
  }
  # secondary x axis
  if (!is.null(topx)) {
    x2 <- data[, topx[1]]
    x1 <- x
    if (-Inf %in% x) x1 <- x[-which(x==-Inf)]
    if (Inf %in% x) x1 <- x[-which(x==Inf)]
    tick <- axisTicks(c(min(x1),max(x1)), log=F)
    x1 <- sort(x1)
    ind <- match(x1[findInterval(tick,x1)], x)
    axis(side=3, at=x[ind], labels=round(x2[ind]))
    mtext(topx[2], side=3, line=0.25, adj=1 )
  }
  # title
  if (!is.null(title)) {
    line <- if (!is.null(topx)) 1.5 else 1
    title(title, line=line)
  }
}


# Box plots 
plot.box <- function(data, ynames, ylab=NULL, 
                     x1name=NULL, x1lab=NULL, x1col=NULL, x1levels=NULL,
                     x2name=NULL, x2lab=NULL, x2col=NULL, x2levels=NULL,
                     znames=NULL, zlab=NULL, zcol=NULL, 
                     hline=NULL, title=NULL, xlab.rot=FALSE, ...) {
  # INPUT:
  # data........data frame where each variable is a column
  # ynames......vector of y variable names, each y is numeric
  #             summary of each y represented by box and whiskers
  #             if x1name=NULL then each y is side by side on the same plot
  #             otherwise separate plot for each y variable
  # ylab........vector of labels for y variables
  #             single value if xnames=NULL
  #             otherwise vector of values corresponding to ynames
  # x1name......x1 is a factor to split by
  #             summary of y at each level of x1
  # x1lab.......label for x1
  # x1col.......vector of color names for each level of x1
  #             ignored if x2name or znames given
  # x1levels....labels for x1 levels
  # x2name......x2 is a second factor to split by (ignored if x1=NULL)
  #             summary of y at each level of x1 split by x2
  # x2lab.......label for x2
  # x2col.......vector of color names for each level of x2
  #             ignored if znames given
  # x2levels....labels for x2 levels
  # znames......z variable name, binary
  #             the proportion of z is calculated over nrow(data)
  #             each box is colored according to the proportion
  #             0 (light) - 1 (dark)
  #             if xnames=NULL, vector length of ynames
  #             otherwise single variable, proportion split by levels of x
  # zlab........label for z variable
  # zcol........color name for z variable
  # hline.......draw a dotted horizontal line
  #             single value if xnames=NULL
  #             otherwise vector of values corresponding to ynames
  # title.......plot title
  # xlab.rot....logical, rotate x labels perpindicular to axis
  # ... other arguments to function boxplot
  #
  # OUTPUT:
  # various box plots

  
  # graphics parameters
  par(mar=c(3,3,1,1)+0.1, mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  if (!is.null(title)) par(mar=par("mar")+c(0,0,1,0))
  
  col <- NULL
  xlab <- NULL
  xnames <- NULL
  legmain <- ""
  if (!is.null(x2name) | !is.null(znames)) x1col <- NULL
  if (is.null(x1name) | !is.null(znames)) x2col <- NULL
  
  # x variables
  if (!is.null(x1name)) {
    data[,x1name] <- as.factor(data[,x1name])
    xcall <- x1name
    if (is.null(x1lab)) x1lab <- gsub("."," ", x1name, fixed=T)
    if(is.null(x1levels)) 
      x1levels <- gsub("."," ", levels(data[,x1name]), fixed=T)
    n1 <- length(x1levels)
    if (!is.null(x1col)) {
      col <- x1col
      legcol <- col
      leglab <- x1levels
      legmain <- x1lab
      xnames <- rep("",n1)
    } else {
      xcomb <- data[, x1name]
      xlab <- x1lab
      xnames <- x1levels
      xseq <- 1:n1
    }
    if (!is.null(x2name)) {
      data[,x2name] <- as.factor(data[,x2name])
      xcall <- paste(x2name, "*", x1name)
      if (is.null(x2lab)) x2lab <- gsub("."," ", x2name, fixed=T)
      if(is.null(x2levels)) x2levels <- levels(data[,x2name])
      n2 <- length(x2levels)
      comb <- n1*n2
      if(!is.null(x2col)) {
        col <- rep(x2col,length=comb)
        legcol <- x2col
        leglab <- x2levels
        legmain <- x2lab
        xseq <- (1:n1-1)*n2+(1+n2)/2
      } else {
        xcomb <- interaction(data[, x2name], data[, x1name])
        xlab <- paste(x1lab, x2lab, sep=":")
        xnames <- levels(xcomb)
        xseq <- 1:comb
      }
    } 
    xnames <- gsub("."," ", xnames, fixed=T)
  } else xlab <- gsub("."," ", ynames, fixed=T)
  
  # z variable
  if (!is.null(znames)) {
    source("color.R")
    if(is.null(x1name))
      zprob <- colMeans(data[,znames]) else 
        zprob <- tapply(data[,znames], xcomb, mean)
    probcol <- color.vec(zprob, zcol)
    col <- probcol$vec
    key <- probcol$key
    legcol <- key$color
    leglab <- key$label
    leglab[seq(1,20,by=2)] <- ""
    if (!is.null(zlab)) legmain <- zlab else
      if (is.null(zlab) & length(znames)==1)
        legmain <- gsub("."," ", znames, fixed=T)        
  }
  
  # PLOTS 
  
  if (!is.null(col)) {
    legmax <- max(strwidth(c(paste0("cc",leglab),legmain),"inches"))
    par(mai=par("mai")+c(0,0,0,legmax))
  }
    
  # only y variables
  if (is.null(x1name)) {
    boxplot(data[,ynames], ylab=ylab, col=col, pch=19, names=xlab, ...)
    if (!is.null(znames)) {
      xpos <- par("usr")[2]*1.01
      ypos <- par("usr")[4]
      if (xlab.rot) ypos <- ypos*1.04
      legend(x=xpos, y=ypos, xpd=T, legend=leglab, pch=15, col=legcol, 
             bty="n", pt.cex=3.5, cex=0.8, y.intersp=1.5, title=legmain, 
             xjust=0, x.intersp=1.3)
    }
    if (!is.null(hline)) abline(h=hline,lty=3)
    if (!is.null(title)) title(title)
  }
  
  if (xlab.rot) {
    xmax <- max(strwidth(xnames,"inches"))-0.5
    par(mai=par("mai")+c(xmax,0,0,0), mgp=par("mgp")+c(0,0.2,0))
    las=2
  } else las=0
  if (xlab=="") par(mar=par("mar")-c(1,0,0,0))

  # y and x variables
  if (!is.null(x1name)) {
    if (is.null(ylab)) ylab <- gsub("."," ", ynames, fixed=T)
    for (i in seq_along(ynames)) {
      formula <- as.formula(paste(ynames[i],xcall,sep="~"))
      boxplot(formula, data, ylab=ylab[i], pch=19,
              xaxt="n", col=col, ...)
      if (!is.null(xnames))
        axis(side=1, at=xseq, labels=xnames, las=las, cex=0.8)
      mtext(xlab, side=1, line=par("mar")[1]-2)
      if (!is.null(col)) {
        xpos <- par("usr")[2]*1.01
        ypos <- par("usr")[4]
        if (xlab.rot) ypos <- ypos*1.04
        legend(x=xpos, y=ypos, xpd=T, 
               legend=leglab, pch=15, col=legcol, bty="n",
               pt.cex=3.5, cex=0.8, y.intersp=1.5, title=legmain, 
               xjust=0, x.intersp=1.3)
        
      }
      if (!is.null(hline)) abline(h=hline[i],lty=3)
      if (!is.null(title)) title(title)
    }
  }
}

# Bar plots 
plot.bar <- function(data, ynames, ylab=NULL, ycol=NULL, 
                     xname=NULL, xlab=NULL, labels=NULL,
                     hline=NULL, hlty=3, hlwd=1,
                     title=NULL, lab.rot=FALSE, 
                     beside=TRUE, overlay=FALSE, ...) {
  # INPUT:
  # data........data frame where each variable is a column
  # ynames......vector of y variable names, each y is numeric
  #             mean of each y represented by bars
  #             if xname is given, mean is calculated within each level of x
  # ylab........labels for y axis
  # ycol........vector of color names for each y variable
  # xname.......x is a factor to split by
  # xlab........label for x axis
  # labels......labels for x axis
  #             if xname=NULL labels are for y variables
  #             otherwise labels are for labels of x variable
  # hline.......draw a horizontal line at value
  # hlty........line type of horizontal line
  # hlwd........line width of horizontal line
  # title.......plot title
  # lab.rot.....logical, rotate axes labels perpindicular to axis
  # beside......logical, if TRUE bars are beside each other, else stacked
  # overlay.....logical, place bars on top of each other
  # ... other arguments to function barplot
  #
  # OUTPUT:
  # various bar plots
  
  
  # graphics parameters
  par(mar=c(3,3,1,1)+0.1, mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  if (!is.null(title)) par(mar=par("mar")+c(0,0,1,0))
  
  
  if (is.null(ylab) & length(ynames)==1) 
    ylab <- gsub("."," ", ynames, fixed=T)
  if (is.null(xname)) {
    if(is.null(labels)) labels <- gsub("."," ", ynames, fixed=T)
    bar <- colMeans(data[,ynames]) 
  } else {
    data[,xname] <- as.factor(data[,xname])
    if (is.null(xlab)) xlab <- gsub("."," ", xname, fixed=T)
    if (is.null(labels)) 
      labels <- gsub("."," ", levels(data[,xname]), fixed=T)
    bar <- NULL
    for (i in seq_along(ynames)) 
      bar <- rbind(bar, tapply(data[,ynames[i]],data[,xname],mean))
  } 

  if (!is.null(ycol)) {
    leglab <- gsub("."," ", ynames, fixed=T)
    legmax <- max(strwidth(paste0("ccc",leglab),"inches"))
    par(mai=par("mai")+c(0,0,0,legmax))              
  }
  
  if (lab.rot) {
    xmax <- max(strwidth(labels,"inches"))-0.4
    par(mai=par("mai")+c(xmax,0.4,0,0), mgp=par("mgp")+c(0.5,0.2,0))
    las=2
  } else las=0
  
  if (xlab=="") par(mar=par("mar")-c(1,0,0,0))
  
  if (overlay) {
    barplot(bar[1,], ylab=ylab, names.arg=labels, 
            las=las, col=ycol[1],...)
    for (i in 2:length(ynames))
      barplot(bar[i,], col=ycol[i], add=T, axes=F, axisnames=F, ...)
  } else {
    barplot(bar, beside=beside, ylab=ylab, names.arg=labels, 
            las=las, col=ycol, ...)    
  }
  mtext(xlab, side=1, line=par("mar")[1]-2)
  if (!is.null(ycol))
    legend(x=par("usr")[2]*1.01, y=par("usr")[4], xpd=T, 
           legend=leglab, pch=15, col=ycol, bty="n",
           pt.cex=3.5, cex=0.8, y.intersp=1.5, xjust=0, x.intersp=1.3)
  if (!is.null(hline)) abline(h=hline,lty=hlty,lwd=hlwd)
  if (!is.null(title)) title(title)
}

# Histogram with normal densities
plot.hist <- function(var, head=NULL, xlab=NULL, ylab=NULL, col,
                      dens=F, dens.col=NULL, dens.lwd=2, dens.lty=1,
                      norm=T, norm.col="black", norm.lwd=2, norm.lty=1, ...) {
  # INPUT:
  # var.........variable to plot
  # head........plot title
  # xlab........label for x axis
  # ylab........label for y axis
  # col.........color to fill histogram, will be 50% lighter
  # dens........add density curve?
  # dens.col....color of density curve
  # dens.lwd....line width of density curve
  # dens.lty....line type of density curve
  # norm........add normal density?
  # norm.col....color of normal curve
  # norm.lwd....line width of normal curve
  # norm.lty....line type of normal curve
  # ... other arguments to function hist
  #
  # OUTPUT:
  # histogram for probability of var with density and normal density curve
  
  source("color.R")
  fill <- color.tint(col, alpha=100, tint=-0.5)
  hist(var, prob=T, main=head, xlab=xlab, ylab=ylab, col=fill, ...)
  if (dens) {
    if (is.null(dens.col)) dens.col <- color.tint(col, tint=0.25)
    lines(density(var), col=dens.col, lwd=dens.lwd, lty=dens.lty)
  }
  if (norm) curve(dnorm(x, mean=mean(var), sd=sd(var)), add=TRUE, 
                  col=norm.col, lwd=norm.lwd, lty=norm.lty)
}

# Scatterplot matrix
plot.corr <- function(x, y, col, ...) {
  # INPUT:
  # x...........x matrix
  # y...........y variable
  # col.........color
  # ... other arguments to function pairs
  #
  # OUTPUT:
  # scatterplot matrix 
  
  par(xaxt="n",yaxt="n")
  source("color.R")
  if(missing(col)) col <- "steelblue"
  
  dat <- data.frame(y=y,x)
  pairs(y~., data=dat, gap=0, col=col,
        text.panel=NULL,
        lower.panel=function(x, y, col) {
          usr <- par("usr"); on.exit(par(usr))
          par(usr=c(0,1,0,1),xaxt="n",yaxt="n")
          r <- round(cor(x,y),2)
          rcol <- if(abs(r)>0.5) "white" else "black"
          col <- color.vec(abs(r),col,200)$vec
          rect(0,0,1,1,col=col)
          text(0.5,0.5,r,cex=2,col=rcol)}, 
        diag.panel=function(x,...) {
          par(new = TRUE)
          tab <- table(x)
          if (length(tab) %in% c(2,3,4)) {
            col <- color.vec(tab/sum(tab), col, 100)$vec
            pie(tab, col=col, labels="")
            } else
              plot.hist(x, col=col)
            },
        upper.panel= function(x, y, col) {
          par(new = TRUE)
          tabx <- table(x)
          taby <- table(y)
          if (length(tabx) %in% c(2,3,4)) {
            col <- color.tint(col, 100)
            boxplot(y~x, col=col, outline=F) 
          } else if (length(taby) %in% c(2,3,4)) {
            col <- color.tint(col, 100)
            boxplot(x~y, col=col, outline=F, horizontal=T) 
          } else {
            col <- color.tint(col, 50)
            plot(x,y,type="p",pch=16,col=col) 
          }
        }, ...)
  p <- ncol(x)+1
  x.coords = par('usr')[1:2]
  x.offset = 0.06 * (x.coords[2] - x.coords[1])  
  xrng =  (x.coords[2] - x.coords[1]) - 2*x.offset
  x.width = xrng/p 
  
  y.coords = par('usr')[3:4]
  y.offset = 0.06 * (y.coords[2] - y.coords[1])
  yrng =  (y.coords[2] - y.coords[1]) - 2*y.offset
  y.width = yrng/p  
  
  # x-axis labels
  text(seq(x.coords[1] + x.offset + 0.5*x.width, 
           x.coords[2] - x.offset - 0.5*x.width,
           length.out=p),y.coords[1]*0.5, c("y",xnames), xpd=TRUE,
       adj=c(0.5,0.5), cex=1.2)
  
  # y-axis labels
  text(x.coords[1]*1.03, seq(y.coords[1] + y.offset + 0.5*y.width, 
                             y.coords[2] - y.offset - 0.5*y.width, 
              length.out=p), rev(c("y",xnames)),
       xpd=TRUE, adj=c(0.5, 0.5), cex=1.2, srt=90)  
}