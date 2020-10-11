
if(FALSE){
  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  
  cmyk2hex <- function(C,M,Y,K) {
    
    C <- C / 100.0
    M <- M / 100.0
    Y <- Y / 100.0
    K <- K / 100.0
    
    n.c <- (C * (1-K) + K)
    n.m <- (M * (1-K) + K)  
    n.y <- (Y * (1-K) + K)
    
    r.col <- ceiling(255 * (1-n.c))
    g.col <- ceiling(255 * (1-n.m))
    b.col <- ceiling(255 * (1-n.y))
    
    return(rgb2hex(r.col, g.col, b.col))
    
  }
  
  cmyk2hex(85,50,0,0) # PD1
  cmyk2hex(35,100,35,10) # PD2
  cmyk2hex(100,100,25,25) # PD3
  cmyk2hex(0,0,0,80) # C1
  cmyk2hex(0,0,0,40) # C2
  
}