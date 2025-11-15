library(terra)

r2 <- rast(xmin=0, xmax=100, ymin=0, ymax=100, res=1, vals=2)
r1 <- rast(xmin=-20, xmax=100, ymin=20, ymax=80, res=1, vals=1)



m <- mosaic(r1,r2)


plot(m)

#perform intersect
er1 <- ext(r1)
er2 <- ext(r2)
#function to determine relationships
#...
overlap <- (er1[1]  < er2[2]|er2[1]  < er1[2])&(er1[3]  < er2[4]|er2[3]  < er1[4])
if(!overlap){
  r <- merge(r1,r2)
  return(r)
}


r2 <- rast(xmin=0, xmax=100, ymin=0, ymax=100, res=1, vals=2)
r1 <- rast(xmin=-20, xmax=100, ymin=20, ymax=80, res=1, vals=1)



m <- mosaic(r1,r2)

plot(m)

#perform intersect
er1 <- ext(r1)
er2 <- ext(r2)
overlap <- (er1[1]  < er2[2]|er2[1]  < er1[2])&(er1[3]  < er2[4]|er2[3]  < er1[4])
if(!overlap){
  r <- merge(r1,r2)
  return(r)
}
ei <- terra::intersect(ext(r1),ext(r2))
r1i <- r1 |> crop(ei)
r2i <- r2 |> crop(ei)

msk0 <- r1i*0+1

#margins of intersect
mwidth = 1/4
bbrk <- ei[3] + (ei[4]-ei[3])*mwidth
tbrk <- ei[3] + (ei[4]-ei[3])*(1-mwidth)
lbrk <- ei[1] + (ei[2]-ei[1])*mwidth
rbrk <- ei[1] + (ei[2]-ei[1])*(1-mwidth)

bflank <- crop(msk0, ext(ei[1],ei[2],ei[3],bbrk))
bcore <- bflank
ycore <- crop(msk0, ext(ei[1],ei[2],bbrk,tbrk))
tflank <- crop(msk0, ext(ei[1],ei[2],tbrk,ei[4]))
tcore <- tflank

lflank <- crop(msk0, ext(ei[1],lbrk,ei[3],ei[4]))
lcore <- lflank
xcore <- crop(msk0, ext(lbrk,rbrk,ei[3],ei[4]))
rflank <- crop(msk0, ext(rbrk,ei[2],ei[3],ei[4]))
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

msk <- min((y.msk+0.01)/(r.msk+0.01),l.msk)
msk <- min((y.msk+0.01)/(l.msk+0.01),r.msk)
msk <- min((x.msk+0.01)/(t.msk+0.01),b.msk)
msk <- min((x.msk+0.01)/(b.msk+0.01),t.msk)

msk <- x.msk
msk <- y.msk
msk <- min(x.msk,y.msk)
plot(msk)

#full intersect
l.full <- msk0
r.full <- msk0
b.full <- msk0
t.full <- msk0
values(l.full) <- rep(seq(1, 0, length.out = ncol(l.full)), times = nrow(l.full))
values(r.full) <- rep(seq(0, 1, length.out = ncol(r.full)), times = nrow(r.full))
values(b.full) <- rep(seq(0, 1, length.out = nrow(b.full)),each = ncol(b.full))
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

plot(msk)
