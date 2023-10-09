
#rearrange a data set that has consecutive monthly columns  into a single column for each of precipitation, mean daily high, low, and mean temperature. Specify prefix for column that represents each property, assuming each column name follow with a two digit month. May either specify High and low temperature or mean and low temperature and the rest of the temperature attributes will be calculated.
climateLong <- function(x, name = 'name',lat = "lat", lon = "lon", elev = "elev", year = 'year',p.select = 'p', t.select = 't',th.select = 'th', tl.select = 'tl'){
  #assemble the header data for a given station
  if(name %in% colnames(x)){name = x[1,name]}
  if(lat %in% colnames(x)){lat = x[1,lat]}
  if(lon %in% colnames(x)){lon = x[1,lon]}
  if(elev %in% colnames(x)){elev = x[1,elev]}
  if(year %in% colnames(x)){year = x[,year]}
  noTs=TRUE#switch to either look at mean temperature or high temperature
  #find out the column index numbers for each parameter
  if(paste0(t.select,"01") %in% colnames(x)){
    t.colrange = grep(paste0("^",t.select,"01$"), colnames(x)):grep(paste0("^",t.select,"12$"), colnames(x))
    noTs=FALSE
  }
  if(paste0(th.select,"01") %in% colnames(x)){
    th.colrange = grep(paste0("^",th.select,"01$"), colnames(x)):grep(paste0("^",th.select,"12$"), colnames(x))
  }
  if(paste0(tl.select,"01") %in% colnames(x)){
    tl.colrange = grep(paste0("^",tl.select,"01$"), colnames(x)):grep(paste0("^",tl.select,"12$"), colnames(x))}
  if(paste0(p.select,"01") %in% colnames(x)){
    p.colrange = grep(paste0("^",p.select,"01$"), colnames(x)):grep(paste0("^",p.select,"12$"), colnames(x))
  }

  header = data.frame(name=name,lat=lat,lon=lon,elev=elev,year=year)
  if(noTs){
    for(i in 1:12){#i=1
      th = x[,th.colrange[i]]
      tl = x[,tl.colrange[i]]
      t = (th+tl)/2
      p = x[,p.colrange[i]]
      mon=i
      ntab0 = cbind(header, mon, t,th, tl,p)
      if(i==1){ntab=ntab0}else{ntab=rbind(ntab,ntab0)}
    }; ntab0=NULL
  }else{
    for(i in 1:12){#i=1
      t = x[,t.colrange[i]]
      tl = x[,tl.colrange[i]]
      th = t*2-tl
      p = x[,p.colrange[i]]
      mon=i
      ntab0 = cbind(header, mon, t,th, tl,p)
      if(i==1){ntab=ntab0}else{ntab=rbind(ntab,ntab0)}
    }; ntab0=NULL
  }
  ntab <- ntab |> arrange(year, mon)
  return(ntab)}
