SoilWaterBudget <- function(p,e,awc){
  S_prev <- 0.5*awc
  s = p*0
  a = p*0
  u = p*0
  d = p*0
  for (j in 1:(floor(12/length(p))+2)){#run cycle 2 or more times to get reasonable initial soil moisture balance
    for (i in seq(1, length(p))) { #i=1
      ## potential water to work with
      s[i] <- min(awc, S_prev + p[i])
      ## actual evapotranspiration "a" with a reduction in available water in proportion to total capacity
      a[i] <- min(s[i], min(1, s[i] / awc) * e[i])
      ## balance applied to soil water content "s" minus surpluses "u"
      s[i] <- S_prev + p[i]- a[i]
      surplus <- max(0, s[i] - awc)
      s[i] <- s[i] - surplus
      u[i] <- surplus
      d[i] <- max(0,awc-s[i])
      S_prev <- s[i]
    }
  }
  return(data.frame(a,s,u,d))}

