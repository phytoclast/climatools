#this generates a 12 month temperature curve from January and July temperatures based on simple sine wave.
generateTemp <- function(t01,t07){
  mon=c(1,2,3,4,5,6,7,8,9,10,11,12)
  t=round((t01-t07)/2*cos((mon-1)/12*2*3.141592)+(t01+t07)/2,1)
  return(t)}

#this generates a generalized precipitation curve from a mean, maximum, minimum,  precipitation, and maximum and minimum months.
generatePpt <- function(p.mean, p.max, p.min, m.max=6, m.min=2){
  mon=c(1,2,3,4,5,6,7,8,9,10,11,12)
  shift1 <- mon-m.min
  shift2 <- mon-m.max
  pshift1 <- (p.max-p.min)/2*cos((shift1+6)/12*2*3.141592)+(p.max+p.min)/2
  pshift2 <- (p.max-p.min)/2*cos((shift2)/12*2*3.141592)+(p.max+p.min)/2
  pwt1 <- ((cos((shift1)/12*2*3.141592)+1)/2)
  pwt2 <- ((cos((shift2)/12*2*3.141592)+1)/2)
  wts1 <- pwt1*(1-pwt2)
  wts2 <- pwt2*(1-pwt1)
  p0 <- ((wts1*pshift1)+(wts2*pshift2))/(wts1+wts2)
  exx<-seq(0.01,5,0.01)
  for(i in 1:length(exx)){#i=1
    dif0 = (mean(p0^exx[i]*max(p0)/max(p0^exx[i])) - p.mean)^2
    if(i==1){
      exdif = exx[1]
      dif=dif0
    }else{
      if(dif>=dif0){
        dif=dif0
        exdif = exx[i]
      }}}
  # rnd <- runif(12, 0.8,1.2)#option to add randomness to rain
  # p0 <- rnd*p0
  p1 <- p0^exdif*max(p0)/max(p0^exdif)
  p <- round(p1*p.mean/mean(p1),1)
  return(p)}
