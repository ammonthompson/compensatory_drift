#' Generate random variate from density function for D.
#'
#' See equations 5 and 6 in Thompson et al. 2016. Genetics
#' @param sig2 is the diffusion rate parameters. Equation 4 Thompson et al. 2016.
random_D <- function(sigmaD2 = 1, D0 = 0, time = 10^3, Dstar = 100, num_sample = 100){
  D_sample <- c()
  for(i in 1:num_sample){
    missed <- TRUE;
    isNeo <- FALSE;
    t <- time;
    diffusionRate <- sigmaD2;
    sigma2Dt <- diffusionRate * t;
    D <- D0;

    #determine if one of the paralogs neofuncitonalized (touched one of the boundaries).
    pNeoL <- pNeo_at_RorL(Dstar, D, diffusionRate, t, "L")
    pNeoR <- pNeo_at_RorL(Dstar, D, diffusionRate, t, "R")
    pNeo <- pNeoR + pNeoL

    if((Dstar-abs(D)) <= 0) isNeo<-TRUE
    if(runif(1) <= pNeo && !isNeo){
      if(runif(1) <= pNeoL/pNeo){
        D <- -Dstar;
      }else{
        D <- Dstar;
      }
      isNeo <- TRUE;
    }

    #Draw a random value for D(t) from N(D,sigmaD2t) if neofuncitonalization is virtually impossible (to speed up calculation)
    #Draw a random value for D(t) if there is a signficant prob of neo, but neo did not occur
    #Accept-reject algorithm: algorithm works by drawing a box around a the distribution left to right and zero to 2*maxprob.
    #Draw random pairs of d and % of 2*maximum from this box to evaluate its relative probability and use that to determine if
    #that value for d hit or missed.
    if(pNeo < 0.0001 && !isNeo){

      D <- D + rnorm(1)*sqrt(sigma2Dt);
      if(D < -Dstar){D <- -Dstar;}
      if(D > Dstar){D <- Dstar;}

    }else if(!isNeo){

      #Set the left and right boundaries as 4 standard deviations from the mean (D) or the absorbing boundaries if greater than.
      if(-4*sqrt(sigma2Dt)+D > -Dstar) left <- D-4*sqrt(sigma2Dt) else left <- -Dstar
      if(4*sqrt(sigma2Dt)+D < Dstar) right <- D+4*sqrt(sigma2Dt) else right <- Dstar

      maxd <- maxFinder(D, sigma2Dt, Dstar, left, right);
      #rejection sampling
      while(missed){
        imagesMaxd<-0;
        imagesd<-0;
        d<-runif(1)*(right-left) + left;
        nMax<-((sqrt(36.8*sigma2Dt)-(d+D))/(4.0*Dstar))+1;
        for(n in 0:(nMax-1)){
          imagesd <- imagesd -
            exp(-(1/(2*sigma2Dt)*(d-(-D+(4*n+2)*Dstar))^2))+
            exp(-(1/(2*sigma2Dt)*(d-(D-(4*n+4)*Dstar))^2))-
            exp(-(1/(2*sigma2Dt)*(d-(-D-(4*n+2)*Dstar))^2))+
            exp(-(1/(2*sigma2Dt)*(d-(D+(4*n+4)*Dstar))^2));
          imagesMaxd <- imagesMaxd -
            exp(-(1/(2*sigma2Dt)*(maxd-(-D+(4*n+2)*Dstar))^2))+
            exp(-(1/(2*sigma2Dt)*(maxd-(D-(4*n+4)*Dstar))^2))-
            exp(-(1/(2*sigma2Dt)*(maxd-(-D-(4*n+2)*Dstar))^2))+
            exp(-(1/(2*sigma2Dt)*(maxd-(D+(4*n+4)*Dstar))^2));
        }

        pdfnearDmax <- (1/sqrt(2.0*pi*sigma2Dt))*(exp(-(1/(2.0*sigma2Dt))*(maxd-D)^2)+imagesMaxd);
        pdfatd <- (1/sqrt(2.0*pi*sigma2Dt))*(exp(-(1/(2.0*sigma2Dt))*(d-D)^2)+imagesd);
        randProb <- runif(1)*pdfnearDmax*2;
        if(pdfnearDmax <= 0 ){
          pdfnearDmax<-1;
          pdfatd<-1
        }
        if(pdfatd <= 0) pdfatd <- pdfnearDmax
        if(randProb < pdfatd){
          missed<-FALSE;
          D <- d;
        }
      }
    }
    D_sample <- c(D_sample, D)
  }
  return( D_sample );

}

maxFinder <- function(D, sigma2t, Dthresh, start, finish){
  pdfatd=0;
  d=start;
  imD = 0;
  nMax = D;
  increment = (finish - start)/100.0;
  for(i in 1:100){
    d =start + i*increment;
    imD=0;
    nMax=((sqrt(36.8*sigma2t)-(d+D))/(4.0*Dthresh))+1;
    for(n in 0:(nMax-1)){
      imD = imD -
        exp(-(1/(2*sigma2t)*(d-(-D+(4*n+2)*Dthresh))^2))+
        exp(-(1/(2*sigma2t)*(d-(D-(4*n+4)*Dthresh))^2))-
        exp(-(1/(2*sigma2t)*(d-(-D-(4*n+2)*Dthresh))^2))+
        exp(-(1/(2*sigma2t)*(d-(D+(4*n+4)*Dthresh))^2));
    }
    if(pdfatd > ((1/sqrt(2.0*pi*sigma2t))*(exp(-(1/(2.0*sigma2t))*(d-D)^2)+imD))){
      break
    }else{
      pdfatd = (1/sqrt(2.0*pi*sigma2t))*(exp(-(1/(2.0*sigma2t))*(d-D)^2)+imD)
    }
  }
  return( d+0.5*increment);
}

pNeo_at_RorL = function(xDstar, xD0, xsigmaD2, xtime, R_L){

  #determine prob of touched Dstar.
  pNeo=0
  neoMaxIndex = (5.0*sqrt(xsigmaD2 * xtime)-xDstar)/(4.0*xDstar)+2;

  if(R_L == "R"){

    for(n in 0:neoMaxIndex){

      pNeo = pNeo + 2*(pnorm(((4*n+3)*xDstar+xD0)/sqrt(xsigmaD2 * xtime)) -
                         pnorm(((4*n+1)*xDstar-xD0)/sqrt(xsigmaD2 * xtime)))
    }

  }else{

    for(n in 0:neoMaxIndex){

      pNeo = pNeo + 2*(pnorm(((4*n+3)*xDstar-xD0)/sqrt(xsigmaD2 * xtime)) -
                         pnorm(((4*n+1)*xDstar+xD0)/sqrt(xsigmaD2 * xtime)))


    }

  }

  return(pNeo)

}

