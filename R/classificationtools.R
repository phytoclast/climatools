#Statistical Classification Tools
#
#
#
#' maxKappa
#'
#' @param actual Vector of actual presence/absences, 0 or 1.
#' @param predicted Vector of predicted presence between 0 and 1.
#'
#' @return Script loops through classification thresholds between 0 and 1, and reports the kappa statistic for the threshold with the most information.
#' @export
#'
#' @examples
maxKappa <- function(actual, predicted){ for(i in 1:99){
  k <- i/100
  Kappa0 <- ModelMetrics::kappa(actual=actual, predicted=predicted, cutoff = k)
  if(i == 1){ maxkappa = k}else{maxkappa = max(Kappa0, maxkappa)}
}
  return(maxkappa)}


#' Find Threshold for Variable separating two classes.
#'
#' @param x Data frame of environmental and classification data.
#' @param class1 Name of column with probabilities for class 1.
#' @param class2 Name of column with probabilities for class 2.
#' @param variable Name of the column with variable of interest.
#'
#' @return Returns a list object identifying the optimal threshold of the variable to separate the two classes.
#' @export
#'
#' @examples
find.threshold <- function(x, class1, class2, variable){
  x <- as.data.frame(x)
  g1 <- x[,class1]
  g2 <- x[,class2]
  v1 <- x[,variable]

  x <- data.frame(v1=v1,g1=g1,g2=g2)
  x <- subset(x, g1 > 0|g2 > 0)
  rm(v1,g1,g2)
  vmax <- max(x$v1)
  vmin <- min(x$v1)
  class1cor <- cor(x$g1, x$v1)
  class2cor <- cor(x$g2, x$v1)
  for(i in 0:100){#i=0
    ip <- i/100
    thr <- ip*(vmax-vmin)+vmin
    x <- x |> mutate(reclass = ifelse(v1 >= thr,1,0))
    sign <- var(x$g1,x$reclass)-var(x$g2,x$reclass)
    xp <- x |> summarise(g1 = sum(g1), g2 = sum(g2), all = g1+g2, p1 = g1/all, p2 = g2/all,
                         gini=1-p1^2-p2^2,
                         ent = -p1*log(ifelse(p1 <= 0, 0.001,p1))-p2*log(ifelse(p2 <= 0, 0.001,p2)))
    xs <- x |> group_by(reclass) |> summarise(g1 = sum(g1), g2 = sum(g2), all = g1+g2, p1 = g1/all, p2 = g2/all,
                                              gini=1-p1^2-p2^2, ent = -p1*log(ifelse(p1 <= 0, 1E-50,p1))-p2*log(ifelse(p2 == 0, 1E-50,p2))) |>
      mutate(gini0 = all/sum(all)*gini, ent0 = all/sum(all)*ent)
    xs <- xs |> summarize(gini = sum(gini0), ent = sum(ent0))

    px <- x |> group_by(reclass) |> summarise(pclass1=sum(g1)/sum(g1+g2),pclass2=sum(g2)/sum(g1+g2))
    p1 <- px$pclass1[ifelse(sign < 0,1,2)]
    p2 <- px$pclass2[ifelse(sign >= 0,1,2)]

    gini_gain = xp$gini - xs$gini
    ent_gain = xp$ent - xs$ent
    if(i==0){
      maxinfo = ent_gain
      finalentropy = xs$ent
      finalgini = xs$gini
      finalsign = sign
      finalthr = thr
      finalp1 = p1
      finalp2 = p2
    }else if(ent_gain >= maxinfo){
      maxinfo=ent_gain
      finalentropy = xs$ent
      finalgini = xs$gini
      finalsign = sign
      finalthr = thr
      finalp1 = p1
      finalp2 = p2}
  }


  result = list(entropy=finalentropy, gini=finalgini, infogain = maxinfo, sign = ifelse(finalsign >= 0, "class1 positive", "class1 negative"), class1 = class1, class2=class2, class1cor = class1cor, class2cor = class2cor, variable = variable, p1=finalp1, p2=finalp2, threshold = finalthr, narrative = paste(
    class1,": ", variable, ifelse(finalsign >= 0, " >= ", " < "),finalthr," | ",
    class2,": ", variable, ifelse(finalsign >= 0, " < ", " >= "),finalthr))

  return(result)}











#' Find Threshold for multiple Variables separating two classes.
#'
#' @param x Data frame of environmental and classification data.
#' @param class1 Name of column with probabilities for class 1.
#' @param class2 Name of column with probabilities for class 2.
#' @param variables Vector of column names with variables of interest.
#'
#' @return Reports a data frame of variable thresholds and variable strength in separating two classes.
#' @export
#'
#' @examples
find.multithreshold <- function(x, class1, class2, variables){
  n = length(variables)
  for(i in 1:n){
    thisvar <- find.threshold(x, class1, class2, variables[i])
    df0=data.frame(
      variable = thisvar$variable,
      threshold = thisvar$threshold,
      class1 = thisvar$class1,
      class2 = thisvar$class2,
      pclass1 = thisvar$p1,
      pclass2 = thisvar$p2,
      sign = thisvar$sign,
      infogain = thisvar$infogain,
      gini = thisvar$gini,
      entropy = thisvar$entropy,
      class1cor = thisvar$class1cor,
      class2cor = thisvar$class2cor,
      narrative = thisvar$narrative)
    if(i==1){df=df0}else{df=rbind(df,df0)}
  }
  return(df)}
