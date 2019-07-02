
package expressiondrift2;

import java.io.File;
import java.io.FileWriter;
import java.util.Random;
import java.io.IOException;
/**
 *
 * @author AmmonThompson
 */
public class Drifter {   
    protected int filenumber;
    protected int mutations;
    protected int mutationsFixed;
    protected double averagemutsize;
    protected double s;
//    public static double N = Double.parseDouble(args[4]);
    protected double W = 1;
    protected int neofunctionalized = 0;
    protected int pseudogenized = 0;
    protected Runtime methodrt = Runtime.getRuntime();
    public Drifter(){
        mutations = 0;
        mutationsFixed = 0;
    } 
    /*
     * diffusionDriftSim makes use of the approximations of in Thompson et al. 2015.
     * Does not have absorbing boundaries. See Class DriftWIthNeofunctionalization.java for 
     * diffusion between boundaries sims.
     */
    public Double[] diffusionDriftSim(Double[] expression, double branchLength, double diffusionRate, double dstar, double theta) throws IOException{
        double A = expression[0].doubleValue() + expression[1].doubleValue();
        double D = expression[0].doubleValue() - expression[1].doubleValue();
        double thetaS = theta;
        double t = branchLength;
        boolean withNeo=true;
        Random randomN = new Random();
        
        if(withNeo){
            D = DriftWithNeofunctionalization.nextDrawD(diffusionRate, D, t, dstar, thetaS);
            A = thetaS;
        }
        else{
            double sigma2Dt = diffusionRate*t;
            A = thetaS;
            D = randomN.nextGaussian()*Math.sqrt(sigma2Dt)+D;
        }
        return new Double[]{new Double((A+D)/2),new Double((A-D)/2)};
    }
    
/*
 * discreteDriftSim method runs simulations with D* boundary. If the simulation passes the threshold neofuncitonalization
 * happened.
 */    
    public Double[] discreteDriftSim(Double[] expression, double branchLength, double simOmega, double simu, double theta, double simN) throws IOException{
        averagemutsize=0;
        double mutsizetotal=0;
        filenumber++;
// ######## Uncomment if you want to get the results of a single simulation through time ######################
        
//        File singleSimFile = new File("C:\\Users\\AmmonThompson\\Desktop\\DATA_RESEARCH\\"
//                + "Expression Study\\Simulations\\Manuscript\\Final_Draft\\Submission_Genetics\\"
//                + "Resubmission\\ABC_redux\\sim_output\\SingleSimFile"+filenumber+".csv");
//        FileWriter singleSimWriter = new FileWriter(singleSimFile);
//        singleSimWriter.append("p1" + "," +"p2" + "," +"t" + "," +"mutsize" + "," +"Ns" +"\n");            
            
        double omega = simOmega;
        double N = simN;
        double u4N = simu*4*N;   // This is the number of expression changing mutations per time unit in a population
        double thetaS = theta;
        Random randomVar = new Random();
        double time = -Math.log(1-randomVar.nextDouble())/u4N;  
        double p1 = expression[0].doubleValue();
        double p2 = expression[1].doubleValue();
        
        W = Math.exp(-Math.pow((p1+p2-thetaS)/omega,2)/2);
        double probFix;
        while(time <= branchLength){                  
            mutations++;
            double p1mutant = p1;
            double p2mutant = p2;
            if(!(p1==0 && p2==0)){
                if((randomVar.nextBoolean() && p1 != 0) || p2 == 0){
                    p1mutant = p1mutant + randomVar.nextGaussian();
                    probFix = getProbFix(p1mutant, p2mutant, W, N, thetaS, omega, 1.0);
                    if(randomVar.nextDouble() <= probFix){
                        mutationsFixed++;
                        mutsizetotal=mutsizetotal+Math.abs(p1-p1mutant);
                        s=Math.exp(-Math.pow((p1mutant+p2mutant-thetaS)/omega,2)/2)/W-1;
//                        singleSimWriter.append(p1 + "," + p2 + "," + time + "," + (p1mutant - p1)+"," + (N*s)+"\n");
                        p1 = p1mutant;
//                        if(Math.abs(p1-p2)>thetaS){p1=0; p2=thetaS; neofunctionalized=1;break;}
                        if(p1 <= 0){p1=0; p2=thetaS; neofunctionalized=1;break;}
                        if(p1 >= thetaS){p1=thetaS; p2=0; neofunctionalized=1;break;}
                    }
                }                        
                else{
                    p2mutant = p2mutant + randomVar.nextGaussian();
                    probFix = getProbFix(p1mutant, p2mutant, W, N, thetaS, omega, 1.0);
                    if(randomVar.nextDouble() <= probFix){
                        mutationsFixed++;
                        mutsizetotal=mutsizetotal+Math.abs(p2-p2mutant);
                        s=Math.exp(-Math.pow((p1mutant+p2mutant-thetaS)/omega,2)/2)/W-1;
//                        singleSimWriter.append(p1 + "," + p2 + "," + time + "," + (p2mutant - p2)+"," + (N*s)+"\n");
                        p2 = p2mutant; 
//                        if(Math.abs(p1-p2)>thetaS){p2=0; p1=thetaS; neofunctionalized=1;break;}
                        if(p2 <= 0){p2=0; p1=thetaS; neofunctionalized=1;break;}
                        if(p2 >= thetaS){p2=thetaS; p1=0; neofunctionalized=1;break;}
                    }
               }
               W = Math.exp(-Math.pow((p1+p2-thetaS)/omega,2)/2);
               if(W==0.0){neofunctionalized=1; break;}
               time = time + (-Math.log(1-randomVar.nextDouble())/u4N);
           }
           else{neofunctionalized=1; break;}
        }
        averagemutsize=mutsizetotal/(double)mutationsFixed;
//        singleSimWriter.close();
        return new Double[]{new Double(p1),new Double(p2)};
    }
    
/*
 * discreteDriftsim2 explicitly simulates neofuncitonalization and pseudogenization as separate mutation events 
 * and does not make use of D* to model loss of a duplicate to ancestral function.
 */    
    public Double[] discreteDriftSim2(Double[] expression, double simOmega, double theta, double simN) throws IOException{
        double omega = simOmega;
        double N = simN;
        double thetaS = theta;
        Random randomVar = new Random();
        double p1 = expression[0].doubleValue();
        double p2 = expression[1].doubleValue();
        double neoFit = 1.001;
        double probNeo = 0.0001;
        double probPseudo = 0.1;
        double randomNum;
        int whichNeoPseudo = 0;
        W = Math.exp(-Math.pow((p1+p2-thetaS)/omega,2)/2);
        double probFix;
        while(true){
            mutations++;
            double p1mutant = p1;
            double p2mutant = p2;
            randomNum = randomVar.nextDouble();
            if(randomNum < probNeo){
                if(randomVar.nextBoolean()){
                    p1mutant=0; whichNeoPseudo = 1;
                }
                else{
                    p2mutant=0; whichNeoPseudo = 2;
                }
                if(randomVar.nextDouble() < getProbFix(p1mutant,p2mutant,W,N,thetaS,omega,neoFit)){
                    neofunctionalized=whichNeoPseudo; System.out.println("neofunctionalized");break;
                }
            }
            else if(randomNum < probPseudo+probNeo){
                if(randomVar.nextBoolean()){
                    p1mutant=0; whichNeoPseudo = 1;
                }
                else{
                    p2mutant=0; whichNeoPseudo = 2;
                }
                if(randomVar.nextDouble() < getProbFix(p1mutant,p2mutant,W,N,thetaS,omega,1.0)){
                    pseudogenized = whichNeoPseudo; System.out.println("pseudo");break;
                }
            }
            else{
                if(randomVar.nextBoolean()){
                    p1mutant = p1mutant + randomVar.nextGaussian();
                    probFix = getProbFix(p1mutant, p2mutant, W, N, thetaS, omega, 1.0);
                    if(randomVar.nextDouble() <= probFix){
                        mutationsFixed++;
                        p1 = p1mutant;
                    }
                }
                else{
                    p2mutant = p2mutant + randomVar.nextGaussian();
                    probFix = getProbFix(p1mutant, p2mutant, W, N, thetaS, omega, 1.0);   
                    if(randomVar.nextDouble() <= probFix){
                        mutationsFixed++;
                        p2 = p2mutant;
                    }
                }
            }    
           W = Math.exp(-Math.pow((p1+p2-thetaS)/omega,2)/2);
//           if(p1<=0){p1=0;System.out.println("pseudo");break;}
//           if(p2<=0){p2=0;System.out.println("pseudo");break;}

        }
        whichNeoPseudo=0;
        return new Double[]{new Double(p1),new Double(p2)};
    }
    public double getAverageMutSize(){return averagemutsize;}
    public void resetMuts(){
        mutations = 0;
        mutationsFixed = 0;
    }
    public Integer[] getMuts(){
        return new Integer[]{new Integer(mutations),new Integer(mutationsFixed)};
    }
    public void resetPopNeo(){
        neofunctionalized = 0;
    }
    public int isPopNeo(){
        return neofunctionalized;
    }
    public void resetPopPsuedo(){pseudogenized=0;}
    public int isPopPseudo(){return pseudogenized;}
    private static double getProbFix(double p1mutant, double p2mutant, double W, double N, double thetaS, double omegasim, double neo){
        double probFix;
        double s = (Math.exp(-Math.pow((p1mutant+p2mutant-thetaS)/omegasim,2)/2)*neo)/W - 1;
        probFix = (1-Math.exp(-2*s))/(1-Math.exp(-4*N*s));
        if(s == 0.0){probFix=1/(2*N);}
        return probFix;
    }

}
