#construct a monthly temperature and precipitation graph with options to show potential evapotranspiration and soil water use, and an option to show a simple bar graph or smoothed curves.

#' Make a monthly climate graph
#'
#' @param t mean monthly temperature (degrees Celsius)
#' @param p mean monthly preciptation (mm)
#' @param th mean daily high temperature (degrees Celsius)
#' @param tl mean daily low temperature (degrees Celsius)
#' @param e potential evapotranspiration (mm) [optional]
#' @param a actual evapotranspiration (mm) [optional]
#' @param u surplus precipitation (mm) [optional]
#' @param t80 80th percentile temperature  [optional](not yet implemented)
#' @param t20 20th percentile temperature  [optional](not yet implemented)
#' @param p80 80th percentile precipitation  [optional](not yet implemented)
#' @param p20 20th percentile precipitation  [optional](not yet implemented)
#' @param name name of the climate station
#' @param bar TRUE/FALSE option to show bar graph (TRUE) or a smoothed curve.
#'
#' @return Climate graph
#' @export
#'
#' @examples selected = subset(climatools::Norms2010, Station_ID %in% 'USW00023234')
#' climtab <- climateLong(selected, name <- 'Station_Name',lat = "Latitude",lon = "Longitude", elev = "Elevation", p.select = 'pp', year='Year_')
#' climtab <- climtab |> group_by(name, lat, lon, elev, mon) |> summarise(across(c("t","th","tl","p"),mean))
#' climtab <- climtab |> mutate(e = GetPET(mon, lat, th, tl, p))
#' budget = SoilWaterBudget(climtab$p, climtab$e, 150)
#' climtab <- climtab |> mutate(s=budget$s,a=budget$a, u=budget$u, d=budget$d)
#'
#' makeClimplot(t=climtab$t,
#' th=climtab$th,
#' tl=climtab$tl,
#' p=climtab$p,
#' e=climtab$e,
#' a=climtab$a,
#' u=climtab$u,
#' name=climtab$name,
#' bar=F)
#'
#' makeClimplot(t=climtab$t,
#' th=climtab$th,
#' tl=climtab$tl,
#' p=climtab$p,
#' e=climtab$e,
#'
#' name=climtab$name,
#' bar=F)
#'
#' makeClimplot(t=climtab$t,
#' th=climtab$th,
#' tl=climtab$tl,
#' p=climtab$p,
#'
#' name=climtab$name,
#' bar=F)
#'
#' makeClimplot(t=climtab$t,
#' th=climtab$th,
#' tl=climtab$tl,
#' p=climtab$p,
#' e=climtab$e,
#'
#' name=climtab$name,
#' bar=T)
#'
makeClimplot = function(t,p,th=NULL,tl=NULL,e=NULL,a=NULL,u=NULL,t80=NULL,t20=NULL,p80=NULL,p20=NULL,name=NULL, bar=FALSE){
  mon=c(1,2,3,4,5,6,7,8,9,10,11,12)
  noE=F;noA=F
  if(is.null(th)){th=t+5}
  if(is.null(tl)){tl=t-5}
  if(is.null(e)){e=0; noE=T}
  if(is.null(a)){a=0; noA=T}
  if(is.null(u)){u=0}
  if(is.null(t80)){t80=0}
  if(is.null(t20)){t20=0}
  if(is.null(p80)){p80=0}
  if(is.null(p20)){p20=0}
  if(is.null(name)){name=NA}

  #assemble as a dataframe
  climtab<- data.frame(mon=mon,
                       t=t,
                       p=p,
                       th=th,
                       tl=tl,
                       e=e,
                       a=a,
                       u=u,
                       t80=t80,
                       t20=t20,
                       p80=p80,
                       p20=p20,
                       name=name)
  climtab12 <- climtab[12,] |> mutate(mon=0)
  climtab01 <- climtab[1,] |> mutate(mon=13)
  climtab <- rbind(climtab12, climtab, climtab01)


    # spline()
    # loess()
    # approx()
    climtabsmooth <- data.frame(mon = seq(0,13,0.1))
    # mod = spline(climtab$t ~ climtab$mon, method = 'natural', n=length(climtab$mon)*5)
    mod = suppressWarnings(loess(climtab$t ~ climtab$mon, span = 1/3))
    climtabsmooth <- climtabsmooth |> mutate(t = predict(mod, newdata = mon),f=ifelse(t<0,t,0))
    mod = suppressWarnings(loess(climtab$p ~ climtab$mon, span = 1/3))
    climtabsmooth <- climtabsmooth |> mutate(p = predict(mod, newdata = mon))
    mod = suppressWarnings(loess(climtab$e ~ climtab$mon, span = 1/3))
    climtabsmooth <- climtabsmooth |> mutate(e = predict(mod, newdata = mon))
    mod = suppressWarnings(loess(climtab$a ~ climtab$mon, span = 1/3))
    climtabsmooth <- climtabsmooth |> mutate(a = predict(mod, newdata = mon), a = ifelse(a>e,e,ifelse(a < 0,0,a)))
    mod = suppressWarnings(loess(climtab$u ~ climtab$mon, span = 1/3))
    climtabsmooth <- climtabsmooth |> mutate(u = predict(mod, newdata = mon), u = ifelse(u>p,p,ifelse(u < 0,0,u)))


  climplot0 <- ggplot(climtab, aes(x=mon)) +
    scale_x_continuous(name= "Month", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),
                       labels=c('01','02','03','04','05','06','07','08','09','10','11','12'))+
    scale_y_continuous(name= "Temperature",breaks=c(-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45),

                       labels=c('-20 (-4)', '-15 (  5)', '-10 (14)', '-5 (23)', '0 (32)', '5 (41)', '10 (50)', '15 (59)', '20 (68)', '25 (77)', '30 (86)', '35 (95)', '40 (104)', '°C (°F)'),
                       sec.axis = sec_axis(trans = ~.*1, name = "Precipitation",breaks=c(0,5,10,15,20,25,30,35,40,45),

                                           labels = c('0', '25   (1)', '50   (2)', '75   (3)', '100 (4)', '125 (5)', '150 (6)', '175 (7)', '200 (8)', 'mm (in)')))+
    theme(legend.position="bottom", legend.text = element_text(size = 5)) +
    scale_fill_manual("Legend", values = c("Temperature" = "yellow",
                                           "Precipitation" = "cyan",
                                           "PET" = "yellow",
                                           "Deficit" = "yellow",
                                           "Frozen" = "white",
                                           "Surplus" = "cyan",
                                           "Soil Recharge" = "blue",
                                           "PET/Recharge" = "purple",
                                           "Actual ET" = "darkcyan",
                                           "Soil Use" = "green")
    )+
    scale_color_manual("",values = c("Temperature" = "black", "Mean" = "red", "Low" = "red", "High"="red","Growth"="darkgreen"))+
    scale_shape_manual("",values = c("Mean" = 19, "Low" = 6, "High"=2))+
    coord_fixed(ratio = 1/9,xlim = c(1,12), ylim = c(-20, 43))+
    labs(title = paste0("Climate of ",name[1]))

  geom_Frz <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=f, fill='Frozen'),alpha = 1,  color="transparent" )
  geom_PET <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=e/5, fill='PET'), alpha = 1,  color="red" )
  geom_TMP <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=t, fill='Temperature'), alpha = 1,  color="red" )
  geom_AET <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=a/5, fill='Soil Use'), alpha = 1,  color="darkgreen" )
  geom_ppt <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=p/5, fill="Surplus"), alpha = 1,  color="blue")
  geom_ppt2 <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=p/5, fill="Precipitation"), alpha = 1,  color="blue")

  geom_rch <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=(p-u)/5, fill='Soil Recharge'), alpha = 1,  color="transparent")
  geom_pptaet <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=pmin(a,p)/5, fill='Actual ET'), alpha = 1,  color="transparent")
  geom_pptpet <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=pmin(e,p)/5, fill='PET/Recharge'), alpha = 1,  color="transparent")
  geom_pptpet2 <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=pmin(e,p)/5), fill='green', alpha = 1,  color="transparent")
  geom_ppttmp <- geom_area(data=climtabsmooth, stat="identity", aes(x=mon, y=pmin(t,p/5)), fill='green', alpha = 1,  color="transparent")

  geom_pptline <- geom_line(data=climtabsmooth, stat="identity",  aes(y=p/5), alpha = 1, color='blue')
  geom_surline <- geom_line(data=climtabsmooth, stat="identity",  aes(y=(p-u)/5), alpha = 1, color='purple')
  geom_PETline <- geom_line(data=climtabsmooth, stat="identity",  aes(y=e/5), alpha = 1, color='red')
  geom_aetline <- geom_line(data=climtabsmooth, stat="identity",  aes(y=a/5), alpha = 1, color='darkgreen')

  geom_temline <- geom_line(data=climtabsmooth, stat="identity",  aes(color= "Temperature", y=t), alpha = 0.5, linewidth=1)
  geom_tpnt <- geom_point(aes(shape='Mean', y=t), color="black")
  geom_tlpnt <- geom_point(aes(shape='Low', y=tl), color="black")
  geom_thpnt <- geom_point(aes(shape='High', y=th), color="black")

  geom_barppt <- geom_bar(data=climtab[2:13,], stat="identity",aes(fill="Precipitation", y=p/5), alpha = 0.85,  color="blue")
  geom_barPET <- geom_bar(data=climtab[2:13,], stat="identity", aes(fill='PET', y=e/5), alpha = 0.60,  color="red" )
  geom_barTMP <- geom_bar(data=climtab[2:13,], stat="identity", aes(fill='PET', y=ifelse(t>0,t,0)), alpha = 0.60,  color="red" )

  geom_perror <- geom_errorbar(aes(ymin=p20/5, ymax=p80/5), width=.2,position=position_dodge(-0.9), color="blue")
  geom_terror <- geom_errorbar(aes(ymin=t20, ymax=t80), width=.2,position=position_dodge(0.9), color="black", alpha=0.5)

  if(bar){
    if(!noE){
      climplot <- climplot0+geom_barppt+geom_barPET+
        geom_temline+geom_tpnt+geom_thpnt+geom_tlpnt
    }else{
      climplot <- climplot0+geom_barppt+geom_barTMP+
        geom_temline+geom_tpnt+geom_thpnt+geom_tlpnt
    }

  }else{
    if(!noA){
      climplot <- climplot0+geom_PET+geom_ppt+geom_rch+geom_AET+geom_pptpet+geom_pptaet+geom_PETline+geom_pptline+geom_aetline+geom_Frz+
        geom_temline+geom_tpnt+geom_thpnt+geom_tlpnt
    }else if(!noE){
      climplot <- climplot0+geom_PET+geom_ppt2+geom_pptpet2+geom_PETline+geom_pptline+geom_Frz+
        geom_temline+geom_tpnt+geom_thpnt+geom_tlpnt
    }else{
      climplot <- climplot0+geom_TMP+geom_ppt2+geom_ppttmp+geom_pptline+geom_Frz+
        geom_temline+geom_tpnt+geom_thpnt+geom_tlpnt
    }}

  return(climplot)
}
