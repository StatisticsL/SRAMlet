wavedoppler <- function(x, snr){
  y = snr * sqrt(x * (1 - x)) * sin((2.1 * pi)/(x + 0.05))/0.289
  return(y)
}

waveblocks <- function(x, snr){
  pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  y = rep(0, length(x));
  for(j in 1:length(pos))
    y = y + ((1 + sign(x - pos[j])) * hgt[j])/2;
  end
  y = (snr * y)/1.914;
  return(y)
}

wavebumps <- function(x, snr){
  pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt = c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  wth = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
  y = rep(0, length(x))
  for(j in 1:length(pos))
    y = y + hgt[j]/(1 + abs(x - pos[j])/wth[j])^4
  end
  y=(snr * y)/0.665
  return(y)
}

waveheavisine <- function(x, snr){
  y = snr * (4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - x))/2.97
  return(y)
}


