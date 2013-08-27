## PhD_Proposal_Spring2013_Figures.r
## Snippets of code for various figures in my PhD research proposal
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: May 15, 2013
## Time-stamp: <2013-08-27 10:33:59 Laura>

setwd("~/Documents/Projects/")
require(colorspace)
#if(!exists("LLbyM")) load("/Users/Laura/Projects/spc-research/SPC_LL_ReadLLData.r")

# returns transparent colours
col2transp <- function(col,tlev=0.5) {

    s1 <- c(col2rgb(col))/255
    s2 <- rgb(s1[1],s1[2],s1[3],alpha=tlev) }

################################################################################################
################################################################################################
## Figures for chapter 2
## Aggregated CPUE for each target species
map.cpue.chap2 <- function(add.map=TRUE) {

    lnow <- LLbyM[LLbyM$flag == "JP",] # japanese fleet only
    lnow$dec <- 10*floor(lnow$yy/10) # assign year to decade

    # pick 6 species with highest overall catch
    top6 <- sort(apply(lnow[,sp.counts],2,sum),decreasing=TRUE)[1:6]
    top6sum <- top6/sum(top6) # proportion of each species in total catch
    top6 <- names(top6)
    top6cpue <- gsub("no","cpue",top6)

    # sum counts + hooks by lat/lon/decade for top 6 species
    d1 <- aggregate(lnow[,c("hhooks",top6)],
                    as.list(lnow[,c("latd","lond","dec")]),
                    sum)
    d1[,top6cpue] <- d1[,top6]/d1$hhooks # calc cpue by 5x5 cell by decade
    d2 <- lapply(top6cpue,function(x) tapply(d1[,x],as.list(d1[,c("lond","latd","dec")]),sum))
    names(d2) <- top6cpue

    # get range of cells by species by decade
    rlond <<- aggregate(d1$lond, as.list(d1[,c("dec","latd")]), range)
    rlatd <- aggregate(d1$latd, as.list(d1[,c("dec","lond")]), range)

    pgrid <- function(wdec) {

            dd.lon <- rlond[rlond$dec == seq(1950,2000,by=10)[wdec],]-2.5
            dd.lat <- rlatd[rlatd$dec == seq(1950,2000,by=10)[wdec],]-2.5

            segments(dd.lon$x.1,dd.lon$latd,
                     dd.lon$x.2,dd.lon$latd)
            segments(dd.lat$lond,dd.lat$x.1,
                     dd.lat$lond,dd.lat$x.2)
    }

    # make maps
    check.dev.size(13.7, 7.65)
    par(family="serif",mai=rep(0.1,4),
        omi=c(0,0,0,0.1))
    layout(matrix(1:49,ncol=7,byrow=TRUE),
           width=c(1,rep(1.5,6)),height=c(0.65,rep(1,6)))

    # barchart at top summarising number of individuals caught that decade:
    sdec <- tapply(apply(d1[,top6],1,sum),list("dec"=d1$dec), sum)

    bcfunk <- function(wdec) {
        cold <- rep("lightyellow2",7)
        cold[wdec] <- "lightskyblue4"
        par(mai=c(0,1,0.2,0.1))
        bh <- barplot(sdec[1:6], names.arg=NA, axes=FALSE, col=cold,border=NA)
        text(bh[wdec],sdec[wdec],round(sdec[wdec]/1000000,0),xpd=NA,pos=1,col="white",
             offset=0.4)

        # add labels
        mtext(seq(1950,2000,by=10)[wdec],line=-0.85,cex=1.5,
              adj=-2.05,side=1,xpd=NA)
    }

    plot.new(); sapply(1:6, bcfunk)
    plotfunk <- function(wsp) {

        dnow <- d2[[wsp]]
        breakv <- c(-9999,
                    seq(0,quantile(dnow,0.95,na.rm=TRUE),length.out=20),
                    9999)
        colv <- rev(heat_hcl(length(breakv)-1))
        dn <- lapply(dimnames(dnow)[1:2],as.numeric)

        subf <- function(wdec) {

        if(wdec==1) { # add pie chart at start
                pp <- which(names(top6sum) == gsub("_cpue","_no",wsp))
                coln <- rep("lightyellow2",length(top6sum))
                coln[pp] <- "lightskyblue4"
                par("mai"=c(0.1,0,0.3,0))
                pie(top6sum, col=coln,labels=NA,border="white")
                par("mai"=rep(0.1,4))

        # species label
        mtext(spnames[gsub("_cpue","",wsp)],xpd=NA,cex=1.15,line=-1)
        }

            image(dn$lond, dn$latd, dnow[,,wdec],axes=FALSE,ann=FALSE,
                  asp=1,breaks=breakv, col=colv)
            #abline(v=dn$lond-2.5,col="white",lwd=0.02)
            #abline(h=dn$latd-2.5,col="white",lwd=0.02)
        pgrid(wdec)
        #if(wsp=="bum_cpue") axis(1)
        #if(wdec==1) axis(2,las=1)

            # add continents
            if(add.map){
                map('world2Hires', regions=unlist(geo), add=T,
            wrap=TRUE,fill=TRUE,col="grey80")
            }
            box()

        }
        sapply(1:6,subf)
    }

    sapply(top6cpue, plotfunk)

    return(d2)
}

################################################################################################
################################################################################################
# 1. Example of area-abundance relationships as a function of home area size
fig2v1 <- function() {
#check.dev.size(14.2, 5.5)
png(file="PhD_ResearchProposal_DiagramAreaAbundanceHomeRange.png",width=14.2/2*480, height=5.5/2*480, res=200)
par(mai=c(0.5,0.5,0.1,0.1),omi=c(0.15,0,0.15,0))
layout(cbind(1,c(2,3),c(4,5),c(6,7)), width=c(100,50,50,50))
col.now="cyan3"
coln2="cyan4"

# First plot: 100 units of a species, randomly dispersed in the environment
vr1 <- runif(100)
plot(1:100, vr1, pch=19, cex=0.75, axes=FALSE, ann=FALSE,col=coln2)
points(1:100, vr1, cex=5, col=col2transp(col.now,0.25), pch=16); box()

# Second plot: 20 units of a species, randomly dispersed, same home range
vr2 <- runif(20)
plot(1:20, vr2, pch=19, cex=0.75, axes=FALSE, ann=FALSE,col=coln2)
points(1:20, vr2, cex=5, pch=16, col=col2transp(col.now,0.25)); box()

# Third plot: 20 units of a species, randomly dispersed, larger home range
plot(1:20, vr2, pch=19, cex=0.75, axes=FALSE, ann=FALSE,col=coln2)
points(1:20, vr2, cex=20,col=col2transp(col.now,0.25),pch=16); box()

# Fourth plot: relationship between area and abundance for second plot
plot(1:10,type="n", pch=19, cex=0.75, xlab="",ylab="",yaxt="n",cex.axis=1.4)
abline(h=5,lwd=3)
mtext("home range size",side=2,line=1)

# Fifth plot: relationship between home range area and abundance for third plot
plot(1:10,type="n", pch=19, cex=0.75, xlab="",ylab="",yaxt="n",cex.axis=1.4)
abline(2,0.75,lwd=3)
mtext("home range size",side=2,line=1)
mtext("time",side=1,line=3)

# Sixth plot: relationship between area and abundance
plot(1:10,type="n", pch=19, cex=0.75, xlab="",ylab="",yaxt="n",cex.axis=1.4)
abline(8,-0.5,lwd=3)
mtext("geographic range size",side=2,line=1)

# Seventh plot: contnd
plot(1:10,type="n", pch=19, cex=0.75, xlab="",ylab="",yaxt="n",cex.axis=1.4)
abline(h=5,lwd=3)
mtext("geographic range size",side=2,line=1)
mtext("time",side=1,line=3)

dev.off()
}

# adds white background to figure before plotting
white.bg <- function() {
    uc <- par("usr")
    polygon(c(uc[1],uc[1],uc[2],uc[2]),
            c(uc[3],uc[4],uc[4],uc[3]), col="white", border=NA, xpd=NA)}

display.coords <- function(center.x=NA,center.y=NA) {

    uu <- par("usr")
    pp <- par("plt")
    ff <- par("fig")

    # coordinates of display device in usr coordinate units
    # get usr coords of current plotting device
    slx <- (uu[2]-uu[1])/(pp[2]-pp[1])
    x0 <- uu[2] - slx*pp[2] # left edge
    x1 <- slx*1+x0 # right edge
    sly <- (uu[4]-uu[3])/(pp[4]-pp[3])
    y0 <- uu[4]-sly*pp[4] # bottom edge
    y1 <- sly*1+y0 # top edge

    # get usr coords of display device from plotting device
    slfx <- (x1-x0)/(ff[2]-ff[1])
    xf0 <- x1 - slfx*ff[2] # left edge
    xf1 <- slfx + xf0 # right edge
    slfy <- (y1-y0)/(ff[4]-ff[3]) # slope y
    yf0 <- y1 - slfy*ff[4] # bottom edge
    yf1 <- slfy + yf0 # top edge

    if(!missing(center.x)) {xf0 <- slfx*ff[1]+xf0
                            xf1 <- slfx*ff[2]+xf0}
    if(!missing(center.y)) {yf0 <- slfy*ff[3]+yf0
                            yf1 <- slfy*ff[4]+yf0}

return(c(xf0,xf1,yf0,yf1))
}


fig2v2 <- function() {

   # check.dev.size(7.5,9)
    dwidth <- 7.5; dheight <- 9
    word.width <- 6.2
    word.height <- word.width * dheight / dwidth
    png(file="PhD_ResearchProposal_DiagramAreaAbundanceHomeRange.png",
        width=word.width/2*480, height=word.height/2*480, res=200)

    par(mfcol=c(3,2),family="serif",mai=c(0.25,0.1,0.25,0.1),
        omi=c(0.1,0.65,0.25,0.65))
    plot.new()

    # add background polygon
    dc <- display.coords(center.y=TRUE)
    corners2poly(dc,xpd=NA,border=NA,col="lightyellow2")

    ################################
    ################################
    # make figures, individuals spread randomly:
    # fig 1
    coln2 <- col2transp("cyan3")
    par(new=TRUE)
    vr1 <- runif(100)
    plot(1:100,vr1,las=1,axes=FALSE,ann=FALSE,
         panel.first=white.bg(),pch=19,cex=7,col=coln2)
    points(1:100, vr1, pch=19, cex=0.75)
    mtext("randomly distributed",adj=0,cex=1.5,line=1)
    box()

    ri.x <- runif(20)
    ri.y <- runif(20)
    plot(ri.x,ri.y,las=1,axes=FALSE,ann=FALSE,
         panel.first=white.bg(),pch=19,cex=7,col=coln2)
    points(ri.x,ri.y,pch=19,cex=0.75)
    box()
    # add trend explanation
    pu <- par("usr")
    plot(ri.x,ri.y,las=1,axes=FALSE,ann=FALSE,
         panel.first=white.bg(),pch=19,cex=20,col=coln2)
    points(ri.x,ri.y,pch=19,cex=0.75)
    box()

    ################################
    ################################
    # figures, individuals aggregated
    rx <- rnorm(100,mean=5,sd=1.5)
    ry <- rnorm(100,mean=5,sd=1.5)
    plot(rx,ry,
         las=1,axes=FALSE,ann=FALSE,
         panel.first=white.bg(),pch=19,
         xlim=c(0,10),ylim=c(0,10),cex=7,col=coln2)
    points(rx,ry,pch=19,cex=0.75)
    box()
    mtext("with aggregation",adj=0,cex=1.5,line=1)
    di.x <- rnorm(20,mean=5,sd=0.5)
    di.y <- rnorm(20,mean=5,sd=0.5)
    plot(di.x,di.y,
         las=1,axes=FALSE,ann=FALSE,
         panel.first=white.bg(),pch=19,cex=7,
         xlim=c(0,10),ylim=c(0,10),col=coln2)
    points(di.x,di.y,pch=19,cex=0.75)
    box()

    # add empty plot because increase home range size does not apply here
    plot.new()
    dev.off()

}

# less annoying wrapper for polygon() function
corners2poly <- function(bottom.x,top.x,bottom.y,top.y,...) {

    vnames <- c("bottom.x","top.x","bottom.y","top.y")

    if(length(bottom.x)==4) {
        allc <- bottom.x
        bottom.x <- allc[1]
        top.x <- allc[2]
        bottom.y <- allc[3]
        top.y <- allc[4] }

    polygon(c(bottom.x,bottom.x,top.x,top.x),
            c(bottom.y,top.y,top.y,bottom.y),
            ...)
}

map.stuff <- function(){
# EEZ land masses to keep for map --
require(xlsx)
data(world2HiresMapEnv) # loads dataset with landlines
ez_table <- read.xlsx("/Users/Laura/Projects/spc-research/NTFSR_DataSummaries_SupportingData.xlsx", 1)
#spc.region <- map('world2Hires', xlim=c(125,260), ylim=c(-45,45),
#                  namesonly=T, plot=F) # subset mapping region
spc.region <- map('world2Hires', namesonly=T, ylim=c(-40,40), plot=F) #
# Keep SPC countries + Indonesia (to show with PNG)
geo <<- lapply(c("New Zealand","USA","USSR","Hawaii",
                "Canada","Australia","Japan",
                "China","Indonesia","Philippines", unique(ez_table$name_map)),
              function(ee) spc.region[grep(ee,spc.region)])
}

c3.RDdiagram <- function() {

    par(mfcol=c(2,2),mai=rep(0.1,4),omi=c(0.65,0.5,0.5,0.1),
        family="serif",tcl=0.35)
    #layout(matrix(1:12,nrow=6),height=c(1,1,3,1,1,3))
    latv <- seq(-3,3,by=0.01)
    # same fishing mortality everywhere,
    # no immigration
   # plot(1:10,rep(1,10),type="l",axes=FALSE);box()
   # plot(1:10,rep(1,10),type="l",axes=FALSE);box()
    plot(latv, dnorm(latv), type="l",ann=FALSE,axes=FALSE,
         lwd=2,col="dark blue",ylim=c(0,0.55));box()
    lines(latv, 0.65*dnorm(latv),col="dodgerblue2",lwd=2)
    axis(2,labels=NA)
    mtext("no immigration",line=2,cex=1.5)
    mtext("even fishing mortality",side=2,line=2,cex=1.5)

    # higher fishing mortality in core,
    # no immigration
#    plot(1:10,rep(1,10),type="l",axes=FALSE);box()
#    plot(1:10,rep(1,10),type="l",axes=FALSE);box()
    plot(latv, dnorm(latv), type="l",axes=FALSE,ann=FALSE,
         col="dark blue",lwd=2,ylim=c(0,0.55));box()
    lines(latv, dnorm(latv)*(0.95-dnorm(latv)),col="dodgerblue2",lwd=2)
    axis(1,at=-3:3,labels=c(NA,"edges",NA,"core",NA,"edges",NA),cex.axis=1.5)
    axis(2,labels=NA)
    mtext("higher fishing mortality in core",side=2,line=2,cex=1.5)

    # same fishing mortality everywhere
    # no immigration
    plot(latv, dnorm(latv), type="l",axes=FALSE,ann=FALSE,
         lwd=2,col="dark blue",ylim=c(0,0.55));box()
    lines(latv, 0.65*dnorm(latv,sd=0.8),col="dodgerblue2",lwd=2)
    mtext("immigration towards core",line=2,cex=1.5)

    # higher fishing mortality in core,
    # immigration towards core
  #  plot(1:10,rep(1,10),type="l",axes=FALSE);box()
  #  plot(1:10,rep(1,10),type="l",axes=FALSE);box()
    plot(latv, dnorm(latv), type="l",axes=FALSE,
         ann=FALSE,lwd=2,col="dark blue",ylim=c(0,0.55));box()
    lines(latv, dnorm(latv,sd=1)*(0.95-dnorm(latv))*(4*dnorm(latv)),col="turquoise",lwd=2)
    lines(latv, dnorm(latv,sd=0.9)*(0.85-dnorm(latv))*2^(1*dnorm(latv,sd=0.75)),col="dodgerblue2",lwd=2)
    axis(1,at=-3:3,labels=c(NA,"edges",NA,"core",NA,"edges",NA),cex.axis=1.5)

}

spatialrezdiagram <- function(wres,draw=TRUE,imp.range,...) {

    ma <- missing(imp.range)

    # create grid definition
    # high res = 0.25 degree
    xv <- seq(0,10,by=wres)
    yv <- seq(0,10,by=wres)

    # create grid
    midp <- expand.grid(xx=xv,yy=yv)
    midp <- midp[midp$xx + wres <=10,]
    midp <- midp[midp$yy + wres <=10,]

    polydraw <- function(bx,by,res,...) {
        polygon(c(bx,bx,bx+res,bx+res),c(by,by+res,by+res,by),...)}

    range.funk <- function(wrange) {
        rangedf <- get(wrange) # get range at time t
        prange <- point.in.polygon(midp$xx, midp$yy, rangedf$x, rangedf$y)

        # making range more irregular if range not imposed by argument 'arange'
        if(ma) {
        ddist <- 0.5 # buffer zone around range for random cells
        add.y <- rep(ddist,length(rangedf$y));add.x <- add.y
        add.y[range0$y <= mean(rangedf$y)] <- -ddist
        add.x[range0$x <= mean(rangedf$x)] <- -ddist
        rrange <- point.in.polygon(midp$xx+ 0.5*wres, midp$yy+ 0.5*wres,
                               rangedf$x+add.x, rangedf$y+add.y)
        rrange[sample(1:length(prange),0.8*length(prange))] <- 0
        arange <- prange+rrange; arange[arange>1] <- 1
    } else { arange <- imp.range[[wrange]] }

        # draw range
        plot(0:10,0:10,type="n",asp=1,xlim=c(0,10+0.5),ylim=c(0,10+0.5),
             ann=FALSE,axes=FALSE)
        colv=c(NA,"turquoise3")
        sapply(1:nrow(midp), function(x)
           polydraw(midp$xx[x],midp$yy[x],
                    col=colv[arange[x]+1],res=wres))
        if(wrange=="range0") mtext(sprintf("%s degree",wres),line=-1,adj=0.1)


        return(arange)
    }

    aa <- NA
    if(draw) {
   # par(mfrow=c(3,3),mai=rep(0.05,4))
    aa <- sapply(c("range0","range1","range2"),range.funk) }

    invisible(list("df"=midp,"inr"=aa))
}

comprez.range <- function(save=FALSE) {

    # 3 x 3 fig, draw high res range first
    # coarse range based on high resolution range
    layout(cbind(matrix(1:9,byrow=TRUE,ncol=3),10:12),
           width=c(3,3,3,1))
    par(mai=c(0.05,0.05,0.1,0.05),omi=c(0,0,0.25,0.2),family="serif")
    r1 <- spatialrezdiagram(wres=0.25) # create high-res range

    # get coarser resolution range
    res2 <- 1
    res3 <- 5
    r2 <- spatialrezdiagram(wres=res2,draw=FALSE)
    r3 <- spatialrezdiagram(wres=res3,draw=FALSE)

    poly2incl <- function(rw,wres,wrange,wdf) {
        # test cell by cell if any of the points in higher resolution
        # are in coarser resolution cell:

        dfnow <- get(wdf)$df
        # dataframe with points in high res range
        phighres <- r1$df[as.logical(r1$inr[,wrange]),]+0.01

        # poly coordinates for
        pxx <- c(dfnow$xx[rw],dfnow$xx[rw],dfnow$xx[rw]+wres,dfnow$xx[rw]+wres)
        pyy <- c(dfnow$yy[rw],dfnow$yy[rw]+wres,dfnow$yy[rw]+wres,dfnow$yy[rw])
        pp <- point.in.polygon(phighres$xx, phighres$yy, pxx, pyy)
        any(pp)}

    rnames <- c("range0","range1","range2")
    pr2 <- lapply(rnames,function(wr)
                  sapply(1:nrow(r2$df),function(x) poly2incl(x,wres=res2,wr,"r2")))
    names(pr2) <- rnames
    spatialrezdiagram(wres=res2, imp.range=pr2)
    pr3 <- lapply(rnames,function(wr)
                  sapply(1:nrow(r3$df),function(x) poly2incl(x,wres=res3,wr,"r3")))
    names(pr3) <- rnames
    spatialrezdiagram(wres=res3, imp.range=pr3)

    # add plot of area change at the end
    sr <- t(sapply(list(as.data.frame(r1$inr),pr2,pr3),function(x) sapply(x,sum)))
    sr <- sr/sr[,1] # calculate prop area left at t=2 and t=3

    pm <- par("mai"); pm[3] <- pm[3]+0.15; pm[1] <- pm[1]+0.08
    par(mai=pm)
    mat <- rbind(sr[1,],apply(sr,2,diff))[,3:1]
    dfunk <- function(dd) {
    bp <- barplot(as.matrix(mat[,dd]),
            col=c("turquoise3","paleturquoise3","grey80"),
            border=NA,axes=FALSE,names.arg=c(""), space=0.5)
    colv <- rep("white",3)
    if(dd < 3) {
    text(rep(bp,3),cumsum(mat[,dd]),pos=1,c("t=3","t=2","t=1"),
         offset=0.3,col=colv,cex=1.15)
 } else { text(bp,1,"t=3",pos=1,col="white",cex=1.15,offset=0.3)}
    abline(h=c(0,1))
if(dd==1) mtext("relative\narea change",cex=0.85,line=0.5)};

    sapply(1:3,dfunk)

    if(save) {
        dev.copy2pdf(file="LTB_ResearchProposal_SpatialResolutionFig.pdf")}
    return(pr2)
}

# arbitrarily drawn range polygon with locator()
range0 <- list(x = c(8.68, 6.56, 3.87, 1.59, 0.41, 0.34, 2.12,
4.79, 9.21, 9.5, 7.89, 6.9, 6.01, 6.66, 8.87, 10.07, 10.26, 10.29,
9.95, 8.97), y = c(1.21, 0.52, 0.78, 1.69, 4.12, 8.05, 9.83,
10.12, 10.12, 8.63, 7.84, 7.4, 6.71, 5.25, 4.86, 4.45, 3.28,
2.37, 1.96, 1.45))

range1 <- list(x = c(6.49, 5.05, 3.71, 2.62, 1.93, 1.74, 1.98,
2.46, 3.15, 5.29, 6.78, 7.45, 6.56, 5.34, 4.76, 4.71, 5.77, 6.85,
7.86, 8.08), y = c(1.38, 1.55, 1.98, 2.85, 4.12, 5.37, 6.66,
7.69, 8.65, 9.01, 9.01, 8.17, 7.16, 6.64, 6.33, 5.63, 5.25, 4.69,
3.97, 2.58))

range2 <- list(x = c(6.57, 5.37, 4.43, 3.39, 2.95, 2.95, 3.23,
3.77, 4.05, 3.83, 4.05, 4.6, 5.09, 5.8, 6.3, 6.79, 6.9, 6.9,
6.79, 6.68), y = c(2.29, 2.67, 3.33, 4.1, 5.08, 6.29, 7.06, 7.6,
7.6, 6.67, 6.01, 5.3, 5.08, 4.86, 4.26, 3.77, 3, 2.78, 2.73,
2.62))
