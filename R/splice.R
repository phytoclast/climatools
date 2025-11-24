
#' Splice Rasters
#'
#' Function blends two terra rasters using smooth gradients. Alternative to mosaic which only averages the area of overlap.
#'
#' @param r1 First raster
#' @param r2 Second raster
#'
#' @returns Raster combination of two rasters.
#' @export
#'
#' @examples rg <- genrast(c(1,-1,-1,1))
#' @examples r1 <- rg[[1]]
#' @examples r2 <- rg[[2]]
#' @examples plot(mosaic(r1,r2))
#' @examples
#' @examples m1 <- spliceraster(r1,r2)
#' @examples plot(m1)
spliceraster <- function(r1,r2){
  require(terra)
  #get extents of input rasters
  er1 <- terra::ext(r1)
  er2 <- terra::ext(r2)
  #determine if overlap and whether to proceed with more complex mosiacing
  overlap <- (er1[1]  < er2[2]|er2[1] < er1[2])&(er1[3]  < er2[4]|er2[3]  < er1[4])
  if(!overlap){
    r <- terra::merge(r1,r2)
    return(r)
  }

  issame <- function(x,e){
    test <- sum(e == x) == 4 | sum(e == x*-1) == 4
    return(test)
  }
  isreverse <- function(x,e){
    test <- sum(e == x*-1) == 4
    return(test)
  }
  e = climatools::rposition(r1,r2)


  #raster completely inside other raster
  xyinside <- (er1[1] >= er2[1] & er1[2] <= er2[2] &
                 er1[3] >= er2[3] & er1[4] <= er2[4])|
    (er1[1] <= er2[1] & er1[2] >= er2[2] &
       er1[3] <= er2[3] & er1[4] >= er2[4])


  #perform intersect
  ei <- terra::intersect(ext(r1),ext(r2))
  r1i <- r1 |> terra::crop(ei)
  r2i <- r2 |> terra::crop(ei)
  r1internal <- !(er1[1] < ei[1] | er1[2] > ei[2] | er1[3] < ei[3] | er1[4] > ei[4])
  r2internal <- !(er2[1] < ei[1] | er2[2] > ei[2] | er2[3] < ei[3] | er2[4] > ei[4])

  #create mask
  msk <- r1i*0+1

  #full intersect routines using entire gradient of overlapping region
  if(!xyinside){

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



    # if(r1l){
    if(issame(c(1,-1,0,0),e)){
      values(msk) <- rep(seq(0, 1, length.out = ncol(msk)), times = nrow(msk))
    }
    # if(r1r){
    if(issame(c(-1,1,0,0),e)){
      values(msk) <- rep(seq(1, 0, length.out = ncol(msk)), times = nrow(msk))
    }
    # if(r1t){
    if(issame(c(0,0,-1,1),e)){
      values(msk) <- rep(seq(1, 0, length.out = nrow(msk)),each = ncol(msk))
    }
    # if(r1b){
    if(issame(c(0,0,1,-1),e)){
      values(msk) <- rep(seq(0, 1, length.out = nrow(msk)),each = ncol(msk))
    }


    # if(r1ne){
    if(issame(c(-1,1,-1,1),e)){
      r.full <- r.full0()
      t.full <- t.full0()
      pmsk1 <- r.full*t.full
      pmsk2 <- (1-r.full)*(1-t.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
      msk <- msk*1.1-0.05
      msk <- ifel(msk > 1,1,ifel(msk < 0,0,msk))
    }
    # if(r1se){
    if(issame(c(-1,1,1,-1),e)){
      r.full <- r.full0()
      b.full <- b.full0()
      pmsk1 <- r.full*b.full
      pmsk2 <- (1-r.full)*(1-b.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
      msk <- msk*1.1-0.05
      msk <- ifel(msk > 1,1,ifel(msk < 0,0,msk))
    }
    # if(r1sw){
    if(issame(c(1,-1,1,-1),e)){
      l.full <- l.full0()
      b.full <- b.full0()
      pmsk1 <- l.full*b.full
      pmsk2 <- (1-l.full)*(1-b.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
      msk <- msk*1.1-0.05
      msk <- ifel(msk > 1,1,ifel(msk < 0,0,msk))
    }
    # if(r1nw){
    if(issame(c(1,-1,-1,1),e)){
      l.full <- l.full0()
      t.full <- t.full0()
      pmsk1 <- l.full*t.full
      pmsk2 <- (1-l.full)*(1-t.full)
      msk <- (pmsk1+0.01)/(pmsk1+pmsk2+0.02)
      msk <- msk*1.1-0.05
      msk <- ifel(msk > 1,1,ifel(msk < 0,0,msk))
    }
    mwidth = 1/2
    bbrk <- ei[3] + (ei[4]-ei[3])*mwidth
    tbrk <- ei[3] + (ei[4]-ei[3])*(1-mwidth)
    lbrk <- ei[1] + (ei[2]-ei[1])*mwidth
    rbrk <- ei[1] + (ei[2]-ei[1])*(1-mwidth)

    bflank <- terra::crop(msk, ext(ei[1],ei[2],ei[3],bbrk))
    bcore <- bflank
    tflank <- terra::crop(msk, ext(ei[1],ei[2],tbrk,ei[4]))
    tcore <- tflank

    lflank <- terra::crop(msk, ext(ei[1],lbrk,ei[3],ei[4]))
    lcore <- lflank
    rflank <- terra::crop(msk, ext(rbrk,ei[2],ei[3],ei[4]))
    rcore <- rflank

    values(lflank) <- rep(seq(0, 1, length.out = ncol(lflank)), times = nrow(lflank))
    values(rflank) <- rep(seq(1, 0, length.out = ncol(rflank)), times = nrow(rflank))
    values(bflank) <- rep(seq(1, 0, length.out = nrow(bflank)), each = ncol(bflank))
    values(tflank) <- rep(seq(0, 1, length.out = nrow(tflank)), each = ncol(tflank))

    values(lcore) <- 1
    values(rcore) <- 1
    values(bcore) <- 1
    values(tcore) <- 1

    y.msk0 <- function(){
      y.msk <- terra::merge(bflank,tflank)
      return(y.msk)}
    x.msk0 <- function(){
      x.msk <- terra::merge(lflank,rflank)
      return(x.msk)}
    r.msk0 <- function(){
      r.msk <- terra::merge(lcore,rflank)
      return(r.msk)}
    l.msk0 <- function(){
      l.msk <- terra::merge(lflank,rcore)
      return(l.msk)}
    t.msk0 <- function(){
      t.msk <- terra::merge(bcore,tflank)
      return(t.msk)}
    b.msk0 <- function(){
      b.msk <- terra::merge(bflank,tcore)
      return(b.msk)}


    #band surpasses on one side
    # if(xbandr){
    if(issame(c(0,1,-1,-1),e)){
      msk1 <- min((b.full0()*t.full0()+0.01)/(l.full0()+0.01),1)
      msk2 <- min((l.full0()+0.01)/(b.full0()*t.full0()+0.01),1)
      msk <- ((msk1 + 1-msk2)/2)
      if(isreverse(c(0,1,-1,-1),e)){msk <- 1-msk}}

    # if(xbandl){
    if(issame(c(1,0,-1,-1),e)){
      msk1 <- min((b.full0()*t.full0()+0.01)/(r.full0()+0.01),1)
      msk2 <- min((r.full0()+0.01)/(b.full0()*t.full0()+0.01),1)
      msk <- ((msk1 + 1-msk2)/2)
      if(isreverse(c(1,0,-1,-1),e)){msk <- 1-msk}}

    # if(ybandt){
    if(issame(c(-1,-1,0,1),e)){
      msk1 <- min((r.full0()*l.full0()+0.01)/(b.full0()+0.01),1)
      msk2 <- min((b.full0()+0.01)/(r.full0()*l.full0()+0.01),1)
      msk <- ((msk1 + 1-msk2)/2)
      if(isreverse(c(-1,-1,0,1),e)){msk <- 1-msk}}

    # if(ybandb){
    if(issame(c(-1,-1,1,0),e)){
      msk1 <- min((r.full0()*l.full+0.01)/(t.full0()+0.01),1)
      msk2 <- min((t.full0()+0.01)/(r.full0()*l.full0()+0.01),1)
      msk <- ((msk1 + 1-msk2)/2)
      if(isreverse(c(-1,-1,1,0),e)){msk <- 1-msk}}



    #incomplete tongue overlap
    # if(xtongr){
    if(issame(c(-1,1,-1,-1),e)){
      msk <- min((b.full0()*t.full0()+0.01)/(l.full0()+0.01),r.full0())
      if(isreverse(c(-1,1,-1,-1),e)){msk <- 1-msk}}

    # if(xtongl){
    if(issame(c(1,-1,-1,-1),e)){
      msk <- min((b.full0()*t.full0()+0.01)/(r.full0()+0.01),l.full0())
      if(isreverse(c(1,-1,-1,-1),e)){msk <- 1-msk}}

    # if(ytongt){
    if(issame(c(-1,-1,-1,1),e)){
      msk <- min((l.full0()*r.full0()+0.01)/(b.full0()+0.01),t.full0())
      if(isreverse(c(-1,-1,-1,1),e)){msk <- 1-msk}}

    # if(ytongb){
    if(issame(c(-1,-1,1,-1),e)){
      msk <- min((l.full0()*r.full0()+0.01)/(t.full0()+0.01),b.full0())
      if(isreverse(c(-1,-1,1,-1),e)){msk <- 1-msk}}


    #cross
    # if(xycross){
    if(issame(c(1,1,-1,-1),e)){
      y.msk <- y.msk0()
      x.msk <- x.msk0()
      msk1 <- min((y.msk+0.01)/(x.msk+0.01),1)
      msk2 <- min((x.msk+0.01)/(y.msk+0.01),1)
      msk <- (msk1 + 1-msk2)/2
      if(isreverse(c(1,1,-1,-1),e)){msk <- 1-msk}}


    # if(inse){
    if(issame(c(-1,0,0,1),e)){
      msk1 <- min((r.full0()+0.01)/(b.full0()+0.01),1)
      msk2 <- min((b.full0()+0.01)/(r.full0()+0.01),1)
      msk <- (msk1 + 1-msk2)/2
      if(isreverse(c(-1,0,0,1),e)){msk <- 1-msk}}
    # if(insw){
    if(issame(c(0,-1,0,1),e)){
      msk1 <- min((l.full0()+0.01)/(b.full0()+0.01),1)
      msk2 <- min((b.full0()+0.01)/(l.full0()+0.01),1)
      msk <- (msk1 + 1-msk2)/2
      if(isreverse(c(0,-1,0,1),e)){msk <- 1-msk}}
    # if(inne){
    if(issame(c(-1,0,1,0),e)){
      msk1 <- min((r.full0()+0.01)/(t.full0()+0.01),1)
      msk2 <- min((t.full0()+0.01)/(r.full0()+0.01),1)
      msk <- (msk1 + 1-msk2)/2
      if(isreverse(c(-1,0,1,0),e)){msk <- 1-msk}}
    # if(innw){
    if(issame(c(0,-1,1,0),e)){
      msk1 <- min((l.full0()+0.01)/(t.full0()+0.01),1)
      msk2 <- min((t.full0()+0.01)/(l.full0()+0.01),1)
      msk <- (msk1 + 1-msk2)/2
      if(isreverse(c(0,-1,1,0),e)){msk <- 1-msk}}


    if(issame(c(0,-1,1,-1),e)){
      msk <- min((l.full0()+0.01)/(t.full0()+0.01),b.full0())
      if(isreverse(c(0,-1,1,-1),e)){msk <- 1-msk}}

    if(issame(c(-1,0,1,-1),e)){
      msk <- min((r.full0()+0.01)/(t.full0()+0.01),b.full0())
      if(isreverse(c(-1,0,1,-1),e)){msk <- 1-msk}}

    if(issame(c(0,-1,-1,1),e)){
      msk <- min((l.full0()+0.01)/(b.full0()+0.01),t.full0())
      if(isreverse(c(0,-1,-1,1),e)){msk <- 1-msk}}

    if(issame(c(-1,0,-1,1),e)){
      msk <- min((r.full0()+0.01)/(b.full0()+0.01),t.full0())
      if(isreverse(c(-1,0,-1,1),e)){msk <- 1-msk}}

    if(issame(c(1,-1,0,-1),e)){
      msk <- min((b.full0()+0.01)/(r.full0()+0.01),l.full0())
      if(isreverse(c(1,-1,0,-1),e)){msk <- 1-msk}}

    if(issame(c(1,-1,-1,0),e)){
      msk <- min((t.full0()+0.01)/(r.full0()+0.01),l.full0())
      if(isreverse(c(1,-1,-1,0),e)){msk <- 1-msk}}

    if(issame(c(-1,1,-1,0),e)){
      msk <- min((t.full0()+0.01)/(l.full0()+0.01),r.full0())
      if(isreverse(c(-1,1,-1,0),e)){msk <- 1-msk}}

    if(issame(c(-1,1,0,-1),e)){
      msk <- min((b.full0()+0.01)/(l.full0()+0.01),r.full0())
      if(isreverse(c(-1,1,0,-1),e)){msk <- 1-msk}}

    if(issame(c(1,-1,0,-1),e)){
      msk <- min((b.full0()+0.01)/(r.full0()+0.01),l.full0())
      if(isreverse(c(1,-1,0,-1),e)){msk <- 1-msk}}



  }else{

    mwidth = 1/4
    bbrk <- ei[3] + (ei[4]-ei[3])*mwidth
    tbrk <- ei[3] + (ei[4]-ei[3])*(1-mwidth)
    lbrk <- ei[1] + (ei[2]-ei[1])*mwidth
    rbrk <- ei[1] + (ei[2]-ei[1])*(1-mwidth)

    bflank <- terra::crop(msk, ext(ei[1],ei[2],ei[3],bbrk))
    bcore <- bflank
    ycore <- terra::crop(msk, ext(ei[1],ei[2],bbrk,tbrk))
    tflank <- terra::crop(msk, ext(ei[1],ei[2],tbrk,ei[4]))
    tcore <- tflank

    lflank <- terra::crop(msk, ext(ei[1],lbrk,ei[3],ei[4]))
    lcore <- lflank
    xcore <- terra::crop(msk, ext(lbrk,rbrk,ei[3],ei[4]))
    rflank <- terra::crop(msk, ext(rbrk,ei[2],ei[3],ei[4]))
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
      y.msk <- terra::merge(bflank,ycore,tflank)
      return(y.msk)}
    x.msk0 <- function(){
      x.msk <- terra::merge(lflank,xcore,rflank)
      return(x.msk)}
    r.msk0 <- function(){
      r.msk <- terra::merge(lcore,xcore,rflank)
      return(r.msk)}
    l.msk0 <- function(){
      l.msk <- terra::merge(lflank,xcore,rcore)
      return(l.msk)}
    t.msk0 <- function(){
      t.msk <- terra::merge(bcore,ycore,tflank)
      return(t.msk)}
    b.msk0 <- function(){
      b.msk <- terra::merge(bflank,ycore,tcore)
      return(b.msk)}



    # if(xyinside){
    if(issame(c(-1,-1,-1,-1),e)){
      msk <- min(x.msk0(),y.msk0())
      if(isreverse(c(-1,-1,-1,-1),e)){msk <- 1-msk}}
    #full band or strip overlap
    # if(yband){
    if(issame(c(0,0,-1,-1),e)){
      msk <- y.msk0()
      if(isreverse(c(0,0,-1,-1),e)){msk <- 1-msk}}
    # if(xband){
    if(issame(c(-1,-1,0,0),e)){
      msk <- x.msk0()
      if(isreverse(c(-1,-1,0,0),e)){msk <- 1-msk}}


    #inner tongue
    # if(xintongr){
    if(issame(c(0,-1,-1,-1),e)){
      msk <- min(y.msk0(),r.msk0())
      if(isreverse(c(0,-1,-1,-1),e)){msk <- 1-msk}}
    # if(xintongl){
    if(issame(c(-1,0,-1,-1),e)){
      msk <- min(y.msk0(),l.msk0())
      if(isreverse(c(-1,0,-1,-1),e)){msk <- 1-msk}}
    # if(yintongt){
    if(issame(c(-1,-1,0,-1),e)){
      msk <- min(t.msk0(),x.msk0())
      if(isreverse(c(-1,-1,0,-1),e)){msk <- 1-msk}}
    # if(yintongb){
    if(issame(c(-1,-1,-1,0),e)){
      msk <- min(b.msk0(),x.msk0())
      if(isreverse(c(-1,-1,-1,0),e)){msk <- 1-msk}}

    #inner corners
    # if(inse){
    if(issame(c(-1,0,0,-1),e)){
      msk <- min(t.msk0(),l.msk0())
      if(isreverse(c(-1,0,0,-1),e)){msk <- 1-msk}}

    # if(insw){
    if(issame(c(0,-1,0,-1),e)){
      msk <- min(t.msk0(),r.msk0())
      if(isreverse(c(0,-1,0,-1),e)){msk <- 1-msk}}

    # if(inne){
    if(issame(c(-1,0,-1,0),e)){
      msk <- min(b.msk0(),l.msk0())
      if(isreverse(c(-1,0,-1,0),e)){msk <- 1-msk}}

    # if(innw){
    if(issame(c(0,-1,-1,0),e)){
      msk <- min(b.msk0(),r.msk0())
      if(isreverse(c(0,-1,-1,0),e)){msk <- 1-msk}}




    #inner edges
    # if(inr){
    if(issame(c(-1,0,0,0),e)){
      msk <- l.msk0()
      if(isreverse(c(-1,0,0,0),e)){msk <- 1-msk}}
    # if(inl){
    if(issame(c(0,-1,0,0),e)){
      msk <- r.msk0()
      if(isreverse(c(0,-1,0,0),e)){msk <- 1-msk}}
    # if(int){
    if(issame(c(0,0,-1,0),e)){
      msk <- b.msk0()
      if(isreverse(c(0,0,-1,0),e)){msk <- 1-msk}}
    # if(inb){
    if(issame(c(0,0,0,-1),e)){
      msk <- t.msk0()
      if(isreverse(c(0,0,0,-1),e)){msk <- 1-msk}}

  }

  gmsk <- r1i*msk+r2i*(1-msk)

  r <- terra::merge(gmsk,r1,r2)
  return(r)
}





#' Raster positions relative to each other
#'
#' The function assesses the relative positions of two overlapping rasters expressed as a 4 value vector.
#'
#' @param r1 First raster
#' @param r2 Second raster
#'
#' @returns Vector of four values corresponding to xmin, xmax, ymin, and ymax of the spatial extent of the area of overlap, wherein a 1 corresponding to where first raster exceeds the extent, while a -1 corresponds to raster 2 exceeding the extent. A 0 values means that both rasters share the same extent value.
#' @export
#'
#' @examples
rposition <-  function(r1,r2){
  er1 <- terra::ext(r1)
  er2 <- terra::ext(r2)
  ei <- terra::intersect(ext(r1),ext(r2))

  l <- ifelse(er1[1] < ei[1],1,ifelse(er2[1] < ei[1],-1,0))
  r <- ifelse(er1[2] > ei[2],1,ifelse(er2[2] > ei[2],-1,0))
  b <- ifelse(er1[3] < ei[3],1,ifelse(er2[3] < ei[3],-1,0))
  t <- ifelse(er1[4] > ei[4],1,ifelse(er2[4] > ei[4],-1,0))
  e <- c(l,r,b,t)
  return(e)
}


#' Generate overlapping dummy rasters.
#'
#' Function is used exclusively for demonstrating splice function examples.
#'
#' @param e 4 value raster showing which raster is to exceed value of overlap extent.
#'
#' @returns List of two rasters that overlap each other.
#' @export
#'
#' @examples
genrast <- function(e) {
  add1 <- ifelse(e >= 0, e,0)
  add2 <- ifelse(e <= 0, -e,0)
  r1 <- terra::rast(xmin=0-add1[1]*10, xmax=50+add1[2]*10, ymin=0-add1[3]*10, ymax=50+add1[4]*10, res=0.2, vals=1, crs='EPSG:4326')
  r2 <- terra::rast(xmin=0-add2[1]*10, xmax=50+add2[2]*10, ymin=0-add2[3]*10, ymax=50+add2[4]*10, res=0.2, vals=2, crs='EPSG:4326')
  return(list(r1,r2))
}
