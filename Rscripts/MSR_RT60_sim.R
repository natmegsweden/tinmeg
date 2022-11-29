
sabien <- function(B, L, H, NRC_low, NRC_high, NRC_ratio){
  
  vol = B*L*H
  xA = L*H
  yA = B*H
  zA = B*L
  
  totA = (xA + yA + zA) * 2
  
  lowA = totA*(1-NRC_ratio)*NRC_low
  highA = totA*NRC_ratio*NRC_high
  
  RT60 = 0.161*(vol/(lowA+highA))
  
  print(paste0('Total room area: ', totA, " m2"))
  print(paste0('Total room volume: ', vol, " m3"))
  print(paste0('RT60: ', round(RT60, 2), " sec"))
  
  return(RT60)

}

plot(sabien(2.85, 3.90, 2.40, 0.1, 0.95, seq(0, 1, by = 0.01)),
     xaxt = "n", yaxt = "n",
     xlab = "Percentage absorbents",
     ylab = "RT60 (sec)")

title("Change of MSR (2.85 x 3.90 x 2.40 m) RT60 \n with surface NRC 0.1 and 0.95")

#axis(1, at = seq(0, 100, by = 10))
axis(1, at = seq(0, 100, by = 10), tck = 1, lty = 2, col = "gray")
axis(2, tck = 1, lty = 2, col = "gray")
box()
