library(terra)
r1 <- rast(xmin=20, xmax=80, ymin=20, ymax=80, res=1, vals=1, crs='EPSG:4326')
r2 <- rast(xmin=0, xmax=80, ymin=20, ymax=80, res=1, vals=2, crs='EPSG:4326')


m <- mosaic(r1,r2)
plot(m)

m1 <- smoothmosaic(r1,r2)
plot(m1)

smoothmosaic <- function(r1,r2){

  #get extents of input rasters
  er1 <- ext(r1)
  er2 <- ext(r2)
  #determine if overlap and whether to proceed with more complex mosiacing
  overlap <- (er1[1]  < er2[2]|er2[1] < er1[2])&(er1[3]  < er2[4]|er2[3]  < er1[4])
  if(!overlap){
    r <- merge(r1,r2)
    return(r)
  }
  
  #determine relationships to feed into appropriate subroutines
  r2Xinside <- er1[1]  < er2[1] & er1[2]  > er2[2]
  r1Xinside <- er2[1]  < er1[1] & er2[2]  > er1[2]
  r2Yinside <- er1[3]  < er2[3] & er1[4]  > er2[4]
  r1Yinside <- er2[3]  < er1[3] & er2[4]  > er1[4]
  anyinside <- r1Yinside | r2Yinside | r1Xinside | r2Xinside
  
  
  tallign  <- er1[4] == er2[4] 
  lallign  <- er1[1] == er2[1] 
  ballign  <- er1[3] == er2[3] 
  rallign  <- er1[2] == er2[2]
  
  yallign  <- tallign & ballign 
  xallign  <- lallign & rallign
  
  #simple one sided relationships
  r1left   <- er1[1] < er2[1] 
  r1right  <- er1[2] > er2[2] 
  r1top    <- er1[4] > er2[4] 
  r1bottom <- er1[3] < er2[3] 
  #corner relationships
  r1ne <- r1right & r1top
  r1se <- r1right & r1bottom
  r1sw <- r1left & r1bottom
  r1nw <- r1left & r1top
  #full side overlap
  r1l  <- r1left & yallign
  r1r  <- r1right & yallign
  r1t  <- r1top & xallign
  r1b  <- r1bottom & xallign
  
  #raster completely inside other raster
  xyinside <- (er1[1] > er2[1] & er1[2] < er2[2] &
                 er1[3] > er2[3] & er1[4] < er2[4])|
    (er1[1] < er2[1] & er1[2] > er2[2] &
       er1[3] < er2[3] & er1[4] > er2[4])
  #full band or strip overlap
  yband <- xallign
  xband <- yallign 
  #band surpasses on one side
  xbandr <- lallign & !yallign  & (er1[2] > er2[2] & r1Yinside| er1[2] < er2[2] & r2Yinside)
  xbandl <- rallign & !yallign  & (er1[1] < er2[1] & r1Yinside| er1[1] > er2[1] & r2Yinside)
  ybandt <- ballign & !xallign  & (er1[4] > er2[4] & r1Xinside| er1[4] < er2[4] & r2Xinside)
  ybandb <- tallign & !xallign  & (er1[3] < er2[3] & r1Xinside| er1[3] > er2[3] & r2Xinside)
  #incomplete tongue overlap
  xtongr <- !yallign  & (er1[2] > er2[2] & r1Yinside| er1[2] < er2[2] & r2Yinside)&
    (er1[1] > er2[1] & r1Yinside| er1[1] < er2[1] & r2Yinside)
  
  xtongl <- !yallign  & (er1[1] < er2[1] & r1Yinside| er1[1] > er2[1] & r2Yinside)&
    (er1[2] < er2[2] & r1Yinside| er1[2] > er2[2] & r2Yinside)
  
  ytongt <- !xallign  & (er1[4] > er2[4] & r1Xinside| er1[4] < er2[4] & r2Xinside)&
    (er1[3] > er2[3] & r1Xinside| er1[3] < er2[3] & r2Xinside)
  
  ytongb <- !xallign  & (er1[3] < er2[3] & r1Xinside| er1[3] > er2[3] & r2Xinside)&
    (er1[4] < er2[4] & r1Xinside| er1[4] > er2[4] & r2Xinside)
  #cross
  xycross <- (er1[1] > er2[1] & er1[2] < er2[2] &
                er1[3] < er2[3] & er1[4] > er2[4])|
    (er1[1] < er2[1] & er1[2] > er2[2] &
       er1[3] > er2[3] & er1[4] < er2[4])
  #innertongue
  xintongr <- !yallign  & lallign  & (er1[2] > er2[2] & r2Yinside| er1[2] < er2[2] & r1Yinside)
  xintongl <- !yallign  & rallign  & (er1[1] < er2[1] & r2Yinside| er1[1] > er2[1] & r1Yinside)
  yintongt <- !xallign  & ballign  & (er1[4] > er2[4] & r2Xinside| er1[4] < er2[4] & r1Xinside)
  yintongb <- !xallign  & tallign  & (er1[3] < er2[3] & r2Xinside| er1[3] > er2[3] & r1Xinside)
  #innercorners
  inse <- !lallign  & rallign & ballign & !tallign
  insw <- lallign  & !rallign & ballign & !tallign
  inne <- !lallign  & rallign & !ballign & tallign
  innw <- lallign  & !rallign & !ballign & tallign
  #inneredges
  inr <- !lallign  & rallign & ballign & tallign
  inl <- lallign  & !rallign & ballign & tallign
  int <- lallign  & rallign & !ballign & tallign
  inb <- lallign  & rallign & ballign & !tallign
  anyinner <- inse|insw|inne|innw|inr|inl|int|inb
  
  #perform intersect
  ei <- terra::intersect(ext(r1),ext(r2))
  r1i <- r1 |> crop(ei)
  r2i <- r2 |> crop(ei)
  
  #create mask
  msk <- r1i*0+1
  
  #full intersect routines using entire gradient of overlapping region
  if(!anyinside & !anyinner){
    
    l.full0 <- function(){
      values(msk) <- rep(seq(1, 0, length.out = ncol(msk)), times = nrow(msk))
      return(msk)
    }
    r.full0 <- function(){
      values(msk) <- rep(seq(0, 1, length.out = ncol(msk)), times = nrow(msk))
      return(msk)
    }
    b.full0 <- function(){
      values(msk) <- rep(seq(0, 1, length.out = nrow(msk)),each = ncol(msk))
      return(msk)
    }
    t.full0 <- function(){
      values(msk) <- rep(seq(1, 0, length.out = nrow(msk)),each = ncol(msk))
      return(msk)
    }
    if(r1l){
      values(msk) <- rep(seq(0, 1, length.out = ncol(msk)), times = nrow(msk))
    }
    if(r1r){
      values(msk) <- rep(seq(1, 0, length.out = ncol(msk)), times = nrow(msk))
    }
    if(r1t){
      values(msk) <- rep(seq(1, 0, length.out = nrow(msk)),each = ncol(msk))
    }
    if(r1b){
      values(msk) <- rep(seq(0, 1, length.out = nrow(msk)),each = ncol(msk))
    }
    
    if(r1ne){
      r.full <- r.full0()
      t.full <- t.full0()
      pmsk1 <- r.full*t.full
      pmsk2 <- (1-r.full)*(1-t.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
    }
    
    if(r1se){
      r.full <- r.full0()
      b.full <- b.full0()
      pmsk1 <- r.full*b.full
      pmsk2 <- (1-r.full)*(1-b.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
    }
    
    if(r1sw){
      l.full <- l.full0()
      b.full <- b.full0()
      pmsk1 <- l.full*b.full
      pmsk2 <- (1-l.full)*(1-b.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
    }
    
    if(r1nw){
      l.full <- l.full0()
      t.full <- t.full0()
      pmsk1 <- l.full*t.full
      pmsk2 <- (1-l.full)*(1-t.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
    }
    
  
    
  }else{
    
    mwidth = 1/4
    bbrk <- ei[3] + (ei[4]-ei[3])*mwidth
    tbrk <- ei[3] + (ei[4]-ei[3])*(1-mwidth)
    lbrk <- ei[1] + (ei[2]-ei[1])*mwidth
    rbrk <- ei[1] + (ei[2]-ei[1])*(1-mwidth)
    
    bflank <- crop(msk, ext(ei[1],ei[2],ei[3],bbrk))
    bcore <- bflank
    ycore <- crop(msk, ext(ei[1],ei[2],bbrk,tbrk))
    tflank <- crop(msk, ext(ei[1],ei[2],tbrk,ei[4]))
    tcore <- tflank
    
    lflank <- crop(msk, ext(ei[1],lbrk,ei[3],ei[4]))
    lcore <- lflank
    xcore <- crop(msk, ext(lbrk,rbrk,ei[3],ei[4]))
    rflank <- crop(msk, ext(rbrk,ei[2],ei[3],ei[4]))
    rcore <- rflank
    
    values(lflank) <- rep(seq(0, 1, length.out = ncol(lflank)), times = nrow(lflank))
    values(xcore) <- 1
    values(rflank) <- rep(seq(1, 0, length.out = ncol(rflank)), times = nrow(rflank))
    values(bflank) <- rep(seq(1, 0, length.out = nrow(bflank)), each = ncol(bflank))
    values(ycore) <- 1
    values(tflank) <- rep(seq(0, 1, length.out = nrow(tflank)), each = ncol(tflank))
    
    values(lcore) <- 1
    values(rcore) <- 1
    values(bcore) <- 1
    values(tcore) <- 1
    
    y.msk0 <- function(){
      y.msk <- merge(bflank,ycore,tflank)
      return(y.msk)}
    x.msk0 <- function(){
      x.msk <- merge(lflank,xcore,rflank)
      return(x.msk)}
    r.msk0 <- function(){
      r.msk <- merge(lcore,xcore,rflank)
      return(r.msk)}
    l.msk0 <- function(){
      l.msk <- merge(lflank,xcore,rcore)
      return(l.msk)}
    t.msk0 <- function(){
      t.msk <- merge(bcore,ycore,tflank)
      return(t.msk)}
    b.msk0 <- function(){
      b.msk <- merge(bflank,ycore,tcore)
      return(b.msk)}
    
    
    
    
    if(xyinside){msk <- min(x.msk0(),y.msk0())}
    #full band or strip overlap
    if(yband){msk <- y.msk0()
    if(r2Yinside){msk <- 1-msk}}
    if(xband){msk <- x.msk0()
    if(r2Xinside){msk <- 1-msk}}
    #band surpasses on one side
    if(xbandr){msk <- min((y.msk0()+0.01)/(r.msk0()+0.01),1)
    if(r2Yinside){msk <- 1-msk}}  
    if(xbandl){msk <- min((y.msk0()+0.01)/(l.msk0()+0.01),1)
    if(r2Yinside){msk <- 1-msk}}
    if(ybandt){msk <- min((x.msk0()+0.01)/(t.msk0()+0.01),1)
    if(r2Xinside){msk <- 1-msk}}
    if(ybandb){msk <- min((x.msk0()+0.01)/(b.msk0()+0.01),1)
    if(r2Xinside){msk <- 1-msk}}
    #incomplete tongue overlap
    if(xtongr){msk <- min((y.msk0()+0.01)/(r.msk0()+0.01),l.msk0())
    if(r2Yinside){msk <- 1-msk}}
    if(xtongl){msk <- min((y.msk0()+0.01)/(l.msk0()+0.01),r.msk0())
    if(r2Yinside){msk <- 1-msk}}
    if(ytongt){msk <- min((x.msk0()+0.01)/(t.msk0()+0.01),b.msk0())
    if(r2Xinside){msk <- 1-msk}}
    if(ytongb){msk <- min((x.msk0()+0.01)/(b.msk0()+0.01),t.msk0())
    if(r2Xinside){msk <- 1-msk}}
    #inner tongue
    if(xintongr){msk <- min(y.msk0(),r.msk0())
    if(r2Yinside){msk <- 1-msk}}
    if(xintongl){msk <- min(y.msk0(),l.msk0())
    if(r2Yinside){msk <- 1-msk}}
    if(yintongt){msk <- min(t.msk0(),x.msk0())
    if(r2Xinside){msk <- 1-msk}}
    if(yintongb){msk <- min(b.msk0(),x.msk0())
    if(r2Xinside){msk <- 1-msk}}
    #inner corners
    if(inse){msk <- min(t.msk0(),l.msk0())
    if(er1[1] < er2[1]){msk <- 1-msk}}
    if(insw){msk <- min(t.msk0(),r.msk0())
    if(er1[2] > er2[2]){msk <- 1-msk}}
    if(inne){msk <- min(b.msk0(),l.msk0())
    if(er1[1] < er2[1]){msk <- 1-msk}}
    if(innw){msk <- min(b.msk0(),r.msk0())
    if(er1[2] > er2[2]){msk <- 1-msk}}
    #inner edges
    if(inr){msk <- l.msk0()
    if(er1[1] < er2[1]){msk <- 1-msk}}
    if(inl){msk <- r.msk0()
    if(er1[2] > er2[2]){msk <- 1-msk}}
    if(int){msk <- b.msk0()
    if(er1[3] < er2[3]){msk <- 1-msk}}
    if(inb){msk <- b.msk0()
    if(er1[4] > er2[4]){msk <- 1-msk}}
    
    
    
    #cross
    if(xycross){
      y.msk <- y.msk0()
      x.msk <- x.msk0()
      if(er1[1] > er2[1]){
        msk <- min((x.msk+0.01)/(y.msk+0.01),1)
      }else{
        msk <- min((y.msk+0.01)/(x.msk+0.01),1)
      }}
    
    
  }
  
  
    gmsk <- r1i*msk+r2i*(1-msk)
  
  r <- merge(gmsk,r1,r2)
  return(r)
}


plot(r)












##########################################


l.full <- msk
values(l.full) <- rep(seq(1, 0, length.out = ncol(l.full)), times = nrow(l.full))

r.full <- msk
values(r.full) <- rep(seq(0, 1, length.out = ncol(r.full)), times = nrow(r.full))

b.full <- msk
values(b.full) <- rep(seq(0, 1, length.out = nrow(b.full)),each = ncol(b.full))

t.full <- msk
values(t.full) <- rep(seq(1, 0, length.out = nrow(t.full)),each = ncol(t.full))






#corner intercept
#ne.corner
y.corner <- t.full
x.corner <- r.full
pmsk1 <- x.corner*y.corner
pmsk2 <- (1-x.corner)*(1-y.corner)
msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
#se.corner
y.corner <- b.full
x.corner <- r.full
pmsk1 <- x.corner*y.corner
pmsk2 <- (1-x.corner)*(1-y.corner)
msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
#sw.corner
y.corner <- b.full
x.corner <- l.full
pmsk1 <- x.corner*y.corner
pmsk2 <- (1-x.corner)*(1-y.corner)
msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
#nw.corner
y.corner <- t.full
x.corner <- l.full
pmsk1 <- x.corner*y.corner
pmsk2 <- (1-x.corner)*(1-y.corner)
msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
msk0 <- msk

#margins of intersect
mwidth = 1/4
bbrk <- ei[3] + (ei[4]-ei[3])*mwidth
tbrk <- ei[3] + (ei[4]-ei[3])*(1-mwidth)
lbrk <- ei[1] + (ei[2]-ei[1])*mwidth
rbrk <- ei[1] + (ei[2]-ei[1])*(1-mwidth)

bflank <- crop(msk, ext(ei[1],ei[2],ei[3],bbrk))
bcore <- bflank
ycore <- crop(msk, ext(ei[1],ei[2],bbrk,tbrk))
tflank <- crop(msk, ext(ei[1],ei[2],tbrk,ei[4]))
tcore <- tflank

lflank <- crop(msk, ext(ei[1],lbrk,ei[3],ei[4]))
lcore <- lflank
xcore <- crop(msk, ext(lbrk,rbrk,ei[3],ei[4]))
rflank <- crop(msk, ext(rbrk,ei[2],ei[3],ei[4]))
rcore <- rflank

values(lflank) <- rep(seq(0, 1, length.out = ncol(lflank)), times = nrow(lflank))
values(xcore) <- 1
values(rflank) <- rep(seq(1, 0, length.out = ncol(rflank)), times = nrow(rflank))
values(bflank) <- rep(seq(1, 0, length.out = nrow(bflank)), each = ncol(bflank))
values(ycore) <- 1
values(tflank) <- rep(seq(0, 1, length.out = nrow(tflank)), each = ncol(tflank))

values(lcore) <- 1
values(rcore) <- 1
values(bcore) <- 1
values(tcore) <- 1

y.msk <- merge(bflank,ycore,tflank)
x.msk <- merge(lflank,xcore,rflank)
r.msk <- merge(lcore,xcore,rflank)
l.msk <- merge(lflank,xcore,rcore)
t.msk <- merge(bcore,ycore,tflank)
b.msk <- merge(bflank,ycore,tcore)

msk <- min((y.msk+0.01)/(r.msk+0.01),1)
msk <- min((y.msk+0.01)/(l.msk+0.01),1)
msk <- min((x.msk+0.01)/(t.msk+0.01),1)
msk <- min((x.msk+0.01)/(b.msk+0.01),1)

#through
msk <- min((y.msk+0.01)/(r.msk+0.01),l.msk)
msk <- min((y.msk+0.01)/(l.msk+0.01),r.msk)
msk <- min((x.msk+0.01)/(t.msk+0.01),b.msk)
msk <- min((x.msk+0.01)/(b.msk+0.01),t.msk)
#inner tongue
msk <- min(y.msk,l.msk)
msk <- min(y.msk,r.msk)
msk <- min(t.msk,x.msk)
msk <- min(b.msk,x.msk)
#inner corners
msk <- min(t.msk,l.msk)
msk <- min(t.msk,r.msk)
msk <- min(b.msk,l.msk)
msk <- min(b.msk,r.msk)

#cross
msk <- min((x.msk+0.01)/(y.msk+0.01),1)
msk <- min((y.msk+0.01)/(x.msk+0.01),1)


msk <- x.msk
msk <- y.msk
msk <- min(x.msk,y.msk)


msk <- min(y.msk,l.msk)
plot(msk)
