
package expressiondrift2;

import java.util.Random;

/**
 *
 * @author AmmonThompson
 */
public class DriftWithNeofunctionalization {
    public static double nextDrawD(double a, double de, double te, double Dstar, double thet){
        boolean missed=true;
        boolean isNeo=false;
        double t = te;
        double diffusionRate = a;
        double sigma2Dt = diffusionRate*t;
        double D=de;
//        double D=-0.5;
        double left, right;
        Random randomN = new Random();
        double d, imagesd, maxd, imagesMaxd, nMax;
        double pdfnearDmax, pdfatd, randProb, pNeo, pNeoR=0, pNeoL=0, neoMaxIndex;
//Set the left and right boundaries as 4 standard deviations from the mean (D) or the absorbing boundaries if greater than.
        if(-4*Math.sqrt(sigma2Dt)+D > -Dstar){left = D-4*Math.sqrt(sigma2Dt);}else{left = -Dstar;}
        if(4*Math.sqrt(sigma2Dt)+D < Dstar){right = D+4*Math.sqrt(sigma2Dt);}else{right = Dstar;}
//determine if one of the paralogs neofuncitonalized (touched one of the boundaries). 
        neoMaxIndex = (5.0*Math.sqrt(sigma2Dt)-Dstar)/(4.0*Dstar)+1;
        for(int n = 0; n<=neoMaxIndex; n++){
            pNeoR = pNeoR + 2*(Phi(((4*n+3)*Dstar+D)/Math.sqrt(sigma2Dt))-Phi(((4*n+1)*Dstar-D)/Math.sqrt(sigma2Dt)));
            pNeoL = pNeoL + 2*(Phi(((4*n+3)*Dstar-D)/Math.sqrt(sigma2Dt))-Phi(((4*n+1)*Dstar+D)/Math.sqrt(sigma2Dt)));
        }
        pNeo = pNeoR+pNeoL;
        if((Dstar-Math.abs(D))<=0.000000001){isNeo=true;}        
        if(randomN.nextDouble()<=pNeo && !isNeo){
            if(randomN.nextDouble()<=pNeoL/pNeo){
                D=-Dstar;
            }
            else{
                D=Dstar;
            }
            isNeo = true;
        }
//Draw a random value for D(t) from N(D,sigmaD2t) if neofuncitonalization is virtually impossible (to speed up calculation)
        if(pNeo < 0.0001 && !isNeo){
            D = D + randomN.nextGaussian()*Math.sqrt(sigma2Dt);
            if(D<-Dstar){D=-Dstar;}
            if(D>Dstar){D=Dstar;}
        }
//Draw a random value for D(t) if there is a signficant prob of neo, but neo did not occur
//Accept-reject algorithm: algorithm works by drawing a box around a the distribution left to right and zero to 2*maxprob. 
//Draw random pairs of d and % of 2*maximum from this box to evaluate its relative probability and use that to determine if 
 //that value for d hit or missed. 
        else if(!isNeo){
            maxd = maxFinder(D, sigma2Dt, Dstar, left, right);
            while(missed){   //rejection sampling
                imagesMaxd=0;
                imagesd=0;
                d=randomN.nextDouble()*(right-left) + left;
                nMax=((Math.sqrt(36.8*sigma2Dt)-(d+D))/(4.0*Dstar))+1;
                for(int n=0; n<nMax; n++){
                    imagesd = imagesd -
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(d-(-D+(4*n+2)*Dstar),2))))+
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(d-(D-(4*n+4)*Dstar),2))))-
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(d-(-D-(4*n+2)*Dstar),2))))+
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(d-(D+(4*n+4)*Dstar),2))));
                    imagesMaxd = imagesMaxd -
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(maxd-(-D+(4*n+2)*Dstar),2))))+
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(maxd-(D-(4*n+4)*Dstar),2))))-
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(maxd-(-D-(4*n+2)*Dstar),2))))+
                            Math.exp(-(1/(2*sigma2Dt)*(Math.pow(maxd-(D+(4*n+4)*Dstar),2))));
                }

                pdfnearDmax = (1/Math.sqrt(2.0*Math.PI*sigma2Dt))*(Math.exp(-(1/(2.0*sigma2Dt))*Math.pow((maxd-D),2))+imagesMaxd);
                pdfatd = (1/Math.sqrt(2.0*Math.PI*sigma2Dt))*(Math.exp(-(1/(2.0*sigma2Dt))*Math.pow((d-D),2))+imagesd);
                randProb = randomN.nextDouble()*pdfnearDmax*2;
                if(pdfnearDmax<=0){pdfnearDmax=1;pdfatd=1;}
                if(pdfatd<=0){pdfatd=pdfnearDmax;}
                if(randProb < pdfatd){
                    missed=false;
                    D = d;
                }
            }
        }
        return D; 
    }
/*
 * Find the mode of the distribution
 */    
    private static double maxFinder(double D, double sigma2t, double Dthresh, double start, double finish){
        //double pdfatd=-100;
        double pdfatd=0;
        double d=start;
        double imD;
        double nMax;
//        int end=(int)((finish-start)/(0.001*sigma2t))+1;
        double increment = (finish - start)/100.0;
        for(int i=0;i<100;i++){
            d =start + i*increment;
            imD=0;
            nMax=((Math.sqrt(36.8*sigma2t)-(d+D))/(4.0*Dthresh))+1;
            for(int n=0;n<nMax;n++){
                imD = imD -
                        Math.exp(-(1/(2*sigma2t)*(Math.pow(d-(-D+(4*n+2)*Dthresh),2))))+
                        Math.exp(-(1/(2*sigma2t)*(Math.pow(d-(D-(4*n+4)*Dthresh),2))))-
                        Math.exp(-(1/(2*sigma2t)*(Math.pow(d-(-D-(4*n+2)*Dthresh),2))))+
                        Math.exp(-(1/(2*sigma2t)*(Math.pow(d-(D+(4*n+4)*Dthresh),2))));
            }
            //if(pdfatd > (1/Math.sqrt(2.0*Math.PI*sigma2t))*(Math.exp(-(1/(2.0*sigma2t))*Math.pow((d-D),2))+imD)){break;}
            //else{pdfatd = (1/Math.sqrt(2.0*Math.PI*sigma2t))*(Math.exp(-(1/(2.0*sigma2t))*Math.pow((d-D),2))+imD);}
            pdfatd = (1/Math.sqrt(2.0*Math.PI*sigma2t))*(Math.exp(-(1/(2.0*sigma2t))*Math.pow((d-D),2))+imD);
        }
        return d+0.5*increment;
    }
    // fractional error in math formula less than 1.2 * 10 ^ -7.
    // although subject to catastrophic cancellation when z in very close to 0
    // from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
    private static double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 + 
                                            t * ( 0.09678418 + 
                                            t * (-0.18628806 + 
                                            t * ( 0.27886807 + 
                                            t * (-1.13520398 + 
                                            t * ( 1.48851587 + 
                                            t * (-0.82215223 + 
                                            t * ( 0.17087277))))))))));
        if (z >= 0){ return  ans;}
        else{        return -ans;}
    }
    // fractional error less than x.xx * 10 ^ -4.
    // Algorithm 26.2.17 in Abromowitz and Stegun, Handbook of Mathematical.
    private static double erf2(double z) {
        double t = 1.0 / (1.0 + 0.47047 * Math.abs(z));
        double poly = t * (0.3480242 + t * (-0.0958798 + t * (0.7478556)));
        double ans = 1.0 - poly * Math.exp(-z*z);
        if (z >= 0){ return  ans;}
        else{        return -ans;}
    }

    // cumulative normal distribution
    // See Gaussia.java for a better way to compute Phi(z)
    private static double Phi(double z) {
        return 0.5 * (1.0 + erf2(z / (Math.sqrt(2.0))));
    }    
}
