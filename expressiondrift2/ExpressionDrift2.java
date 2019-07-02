
package expressiondrift2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;

/**
 *
 * @author AmmonThompson
 */
public class ExpressionDrift2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException{
        ArrayList<Integer> parentNode = new ArrayList();
        ArrayList<Integer> childNode = new ArrayList();
        ArrayList<Double> brLength = new ArrayList();
//        Scanner scan1 = new Scanner(new File("treeDat2.csv"));
//        Scanner scan1 = new Scanner(new File("treeDat.noE.csv"));        
        Scanner scan1 = new Scanner(new File("oneLineage.csv"));
//        Scanner scan1 = new Scanner(new File("treeDat_test.csv"));
//        Scanner scan1 = new Scanner(new File("treeDatnoE_test.csv"));
 
        
        scan1.nextLine();
        Scanner scan2;
        
//        File rawDatFile = new File("C:\\Users\\AmmonThompson\\Desktop\\DATA_RESEARCH\\"
//                + "Expression Study\\Simulations\\Manuscript\\Final_Draft\\Submission_Genetics\\"
//                + "Resubmission\\ABC_redux\\sim_output\\further_inference\\RawDatFile.csv");
//        File rawDatFile = new File("RawDatFile"+args[0]+".csv");
        File rawDatFile = new File("Rawdatest.csv");
        
        try{
            boolean modelIsDiffusion = true;
            boolean sameNode;
            double omega;
            double N;
            double u;
            double theta;
            double pstar;
            double dstar=0.0;
            double gamma;
            double variance_mutation;
            double diffusionRate;
            Random priorRandomSample = new Random();
            Drifter driftObject = new Drifter();
            Double[] startExpression;            
            variance_mutation = 1; 
            
//##############################################################################                
//######### Parameter defaults #####################################################################            

            N = Math.pow(10,2);
            u = Math.pow(10,-3);
            omega = 2;
            theta = 100;
            pstar = 0;
            dstar=theta-2*pstar;
            diffusionRate = 1.543*(Math.pow(omega,3.0)*u)/(Math.pow(N, 3.0/2.0));

//##############################################################################                
//### Input paramaters from the command line. Uncomment to use. #######                
//##############################################################################                
//            N = Double.parseDouble(args[4]);
//            u = Double.parseDouble(args[2]);
//            omega = Double.parseDouble(args[3]);
//            theta = Double.parseDouble(args[5]);
            
//##############################################################################            
            
// Parse through a tree file output csv from the R package ape to get a list of branch lengths for the tree         
            while(scan1.hasNextLine()){
                scan2 = new Scanner(scan1.nextLine());
                scan2.useDelimiter(",");
                parentNode.add(scan2.nextInt());
                childNode.add(scan2.nextInt());
//                brLength.add(scan2.nextDouble()*Math.pow(10,6)*Double.parseDouble(args[6]));   //branch lengths are in units of my
                brLength.add(scan2.nextDouble()*Math.pow(10,6));   //branch lengths are in units of years

            }
            
            HashMap<Integer, Double[]> nodeExpression = new HashMap();
            ArrayList<Double[]> tipExpression = new ArrayList();
            ArrayList<Integer> tipNeo = new ArrayList();
            
            FileWriter rawDatWriter = new FileWriter(rawDatFile);
            
//##############################################################################################################################            
//################################################ Experiment File Outputs ###################################################################            

            
//            rawDatWriter.append("DiffusionRate" +","+"threshold"+","+"sp1"+","+"sp2"+","+"sp3"+","+"sp4"+","+"sp5"+","+"sp6"+","+"sp7"+","+"sp8"+","+"sp9"+","+"sp10"+
//                    ","+"sp1Neo"+","+"sp2Neo"+","+"sp3Neo"+","+"sp4Neo"+","+"sp5Neo"+","+"sp6Neo"+","+"sp7Neo"+","+"sp8Neo"+","+"sp9Neo"+","+"sp10Neo\n");
            
//            rawDatWriter.append("u" + "," + "ThetaS" + "," + "omega" + ","  + "N" + "," + "isNeo"+","+"Mutations"+ "," + 
//                    "Fixed_mutations" +"," +"average_mut_size"+","+"sp1"+","+"sp2"+","+"sp3"+","+"sp4"+","+"sp5"+","+"sp6"+","+"sp7"+","+"sp8"+","+"sp9"+","+"sp10\n");
            
//            rawDatWriter.append("u" + "," + "ThetaS" + "," + "omega" + ","  + "N" + "," + "isNeo"+","+"Mutations"+ "," + 
//                    "Fixed_mutations" +"," +"average_mut_size"+","+"sp1"+","+"sp2"+","+"sp3"+","+"sp4"+","+"sp5"+","+"sp6"+","+"sp7"+","+"sp8\n");

//            rawDatWriter.append("u"+","+"threshold"+ ","+"omega"+","+"sigmam"+ ","+"neo"+","+"pseudo"+","+"Mutations"+","+
//                    "Fixed_mutations"+","+"p1"+","+"p2"+","+"Ne: "+ N + "," + "generations: "+ brLength.get(0) + "\n");
            
            rawDatWriter.append("u"+","+"theta"+ ","+"omega"+","+ "N" + "," + "neo"+","+"mutations"+ "," + "fixed"+","+"averageSize"+","+"p1"+","+"p2"+"\n");

// for ABC redux ################
//            rawDatWriter.append("DiffusionRate" +","+"pstar"+","+"theta"+","+"sp1"+","+"sp2"+","+"sp3"+","+"sp4"+","+"sp5"+","+"sp6"+","+"sp7"+","+"sp8"+","+"sp9"+","+"sp10"+
//                    ","+"sp1Neo"+","+"sp2Neo"+","+"sp3Neo"+","+"sp4Neo"+","+"sp5Neo"+","+"sp6Neo"+","+"sp7Neo"+","+"sp8Neo"+","+"sp9Neo"+","+"sp10Neo\n");
            
//            rawDatWriter.append("DiffusionRate" +","+"pstar"+","+"theta"+","+"sp1"+","+"sp2"+","+"sp3"+","+"sp4"+","+"sp5"+","+"sp6"+","+"sp7"+","+"sp8\n");
            
            
            

//            rawDatWriter.append("u"+","+"Dstar"+ ","+"omega"+","+ "N" + "," + "neo"+","+"mutations"+ "," + "fixed"+","+"averageSize"+","+"p1"+","+"p2"+"\n");
            
            
//############################################################################################################################################################
            long start=System.currentTimeMillis();
            int samplesize=1;
            for(int f=0; f<samplesize; f++){
            int maxRep = 1000;//Integer.parseInt(args[1]);
            
// replicates of simulations on a tree
            for(int rep = 0; rep < maxRep; rep++){
                
// ABC priors
//                diffusionRate = Math.exp(25.328*priorRandomSample.nextDouble())*Math.pow(10,-11);                
//                gamma=Math.exp(9.4*priorRandomSample.nextDouble())*Math.pow(10,-5);
//                theta=1996.0*priorRandomSample.nextDouble()+4.0;
//                pstar=gamma*theta;
//                dstar=theta-2.0*pstar;
                
// Simulate evolution on each branch of the tree starting at the root and progressing to the tips              
                
                startExpression = new Double[]{new Double(theta/2.0), new Double(theta/2.0)};
                nodeExpression.put(parentNode.get(0),startExpression);
                Double[] terminalExpression;
                for(int i = 0; i < parentNode.size(); i++){
                    if(modelIsDiffusion){
                        if((dstar-(Math.abs(nodeExpression.get(parentNode.get(i))[0]-nodeExpression.get(parentNode.get(i))[1]))) > 0.00000001){
                            terminalExpression = driftObject.diffusionDriftSim(nodeExpression.get(parentNode.get(i)),
                                    brLength.get(i), diffusionRate, dstar, theta);
                        }
                        else{
                            terminalExpression=nodeExpression.get(parentNode.get(i));
                        }
                    }
                    else{                    
                        terminalExpression = driftObject.discreteDriftSim(nodeExpression.get(parentNode.get(i)), brLength.get(i), omega, u, theta, N);
//                        terminalExpression = driftObject.discreteDriftSim2(nodeExpression.get(parentNode.get(i)), omega, theta, N);
                        
                    }
                    nodeExpression.put(childNode.get(i), terminalExpression);
                }
// Record tip expression results in the ArrayList tipExpression                
                for(int i = 0; i < childNode.size(); i++){
                    sameNode = false;
                    for(int j = 0; j < parentNode.size(); j++){
                        if(childNode.get(i).equals(parentNode.get(j))){
                            sameNode = true;
                            break;
                        }
                    }
                    if(!sameNode){
                        tipExpression.add(nodeExpression.get(childNode.get(i)));
                    }
                }
// Loop through simulation on tree tip results and store which if any paralog neofunctionalize, ie. was absorbed by boundary                
                for(int i=0;i<tipExpression.size()-1;i++){
                    if((tipExpression.get(i)[0]-tipExpression.get(i)[1])-dstar>=-0.000000001){
                        tipNeo.add(new Integer(1));
                    }
                    else if((tipExpression.get(i)[0]-tipExpression.get(i)[1])+dstar<=0.000000001){
                        tipNeo.add(new Integer(-1));
                    }
                    else{
                        tipNeo.add(new Integer(0));
                    }
                }
//########################################################################################################################################################################                
//                if(((tipNeo.get(4)==1 && tipNeo.get(8)==1)||(tipNeo.get(4)==-1 && tipNeo.get(8)==-1)) && tipNeo.get(0)==0 && tipNeo.get(1)==0 && tipNeo.get(2)==0   && tipNeo.get(5)==0 &&
//                        tipNeo.get(3)==0 && tipNeo.get(6)==0 && tipNeo.get(7)==0 && tipNeo.get(9)==0){
//                rawDatWriter.append(diffusionRate + ","  + threshold + "," +
//                        (tipExpressionReplicates.get(0)[0]-tipExpressionReplicates.get(0)[1])+","+
//                        (tipExpressionReplicates.get(1)[0]-tipExpressionReplicates.get(1)[1])+","+
//                        (tipExpressionReplicates.get(2)[0]-tipExpressionReplicates.get(2)[1])+","+
//                        (tipExpressionReplicates.get(3)[0]-tipExpressionReplicates.get(3)[1])+","+
//                        (tipExpressionReplicates.get(4)[0]-tipExpressionReplicates.get(4)[1])+","+
//                        (tipExpressionReplicates.get(5)[0]-tipExpressionReplicates.get(5)[1])+","+
//                        (tipExpressionReplicates.get(6)[0]-tipExpressionReplicates.get(6)[1])+","+
//                        (tipExpressionReplicates.get(7)[0]-tipExpressionReplicates.get(7)[1])+","+
//                        (tipExpressionReplicates.get(8)[0]-tipExpressionReplicates.get(8)[1])+","+
//                        (tipExpressionReplicates.get(9)[0]-tipExpressionReplicates.get(9)[1])+","+
//                        tipNeo.get(0).intValue()+","+tipNeo.get(1).intValue()+","+tipNeo.get(2).intValue()+","+
//                        tipNeo.get(3).intValue()+","+tipNeo.get(4).intValue()+","+tipNeo.get(5).intValue()+","+
//                        tipNeo.get(6).intValue()+","+tipNeo.get(7).intValue()+","+tipNeo.get(8).intValue()+","+
//                        tipNeo.get(9).intValue()+"\n");
//                
//                        System.out.println("rep number: "+ rep + "    lamTime: " + (System.currentTimeMillis()-start)); start = System.currentTimeMillis();
//
//                }
               
//########################################################################################################################################################################                
//                rawDatWriter.append(u + "," + theta + "," + omega + ","  + N + "," + driftObject.isPopNeo() + "," + driftObject.getMuts()[0].intValue() + ","+ 
//                          driftObject.getMuts()[1].intValue() + ","+driftObject.getAverageMutSize()+ "," +
//                        (tipExpression.get(0)[0]-tipExpression.get(0)[1])+","+
//                        (tipExpression.get(1)[0]-tipExpression.get(1)[1])+","+
//                        (tipExpression.get(2)[0]-tipExpression.get(2)[1])+","+
//                        (tipExpression.get(3)[0]-tipExpression.get(3)[1])+","+
//                        (tipExpression.get(4)[0]-tipExpression.get(4)[1])+","+
//                        (tipExpression.get(5)[0]-tipExpression.get(5)[1])+","+
//                        (tipExpression.get(6)[0]-tipExpression.get(6)[1])+","+
//                        (tipExpression.get(7)[0]-tipExpression.get(7)[1])+","+
//                        (tipExpression.get(8)[0]-tipExpression.get(8)[1])+","+
//                        (tipExpression.get(9)[0]-tipExpression.get(9)[1])+"\n");
//
//                rawDatWriter.append(u + "," + theta + "," + omega + ","  + N + "," + driftObject.isPopNeo() + "," + driftObject.getMuts()[0].intValue() + ","+ 
//                          driftObject.getMuts()[1].intValue() + ","+driftObject.getAverageMutSize()+ "," +
//                        (tipExpression.get(0)[0]-tipExpression.get(0)[1])+","+
//                        (tipExpression.get(1)[0]-tipExpression.get(1)[1])+","+
//                        (tipExpression.get(2)[0]-tipExpression.get(2)[1])+","+
//                        (tipExpression.get(3)[0]-tipExpression.get(3)[1])+","+
//                        (tipExpression.get(4)[0]-tipExpression.get(4)[1])+","+
//                        (tipExpression.get(5)[0]-tipExpression.get(5)[1])+","+
//                        (tipExpression.get(6)[0]-tipExpression.get(6)[1])+","+
//                        (tipExpression.get(7)[0]-tipExpression.get(7)[1])+"\n");
//
                  rawDatWriter.append(u + "," + theta + "," + omega + ","  + N + "," + driftObject.isPopNeo() + "," + driftObject.getMuts()[0].intValue() + ","+ 
                          driftObject.getMuts()[1].intValue() + ","+driftObject.getAverageMutSize()+","+tipExpression.get(0)[0] +","+ tipExpression.get(0)[1]+"\n");
//
//                  rawDatWriter.append(u + "," + theta + "," + omega + ","  + N + "," + driftObject.isPopNeo() + "," + driftObject.isPopPseudo() + "," + driftObject.getMuts()[0].intValue() + ","+ 
//                          driftObject.getMuts()[1].intValue() + ","+tipExpression.get(0)[0] +","+ tipExpression.get(0)[1]+"\n");
  
// for ABC redux ##################
//                if(((tipNeo.get(4)==1 && tipNeo.get(8)==1)||(tipNeo.get(4)==-1 && tipNeo.get(8)==-1)) && tipNeo.get(0)==0 && tipNeo.get(1)==0 && tipNeo.get(2)==0   && tipNeo.get(5)==0 &&
//                        tipNeo.get(3)==0 && tipNeo.get(6)==0 && tipNeo.get(7)==0 && tipNeo.get(9)==0){
//                rawDatWriter.append(diffusionRate + ","  + pstar + "," + theta + ","+
//                        (tipExpression.get(0)[0]-tipExpression.get(0)[1])+","+
//                        (tipExpression.get(1)[0]-tipExpression.get(1)[1])+","+
//                        (tipExpression.get(2)[0]-tipExpression.get(2)[1])+","+
//                        (tipExpression.get(3)[0]-tipExpression.get(3)[1])+","+
//                        (tipExpression.get(4)[0]-tipExpression.get(4)[1])+","+
//                        (tipExpression.get(5)[0]-tipExpression.get(5)[1])+","+
//                        (tipExpression.get(6)[0]-tipExpression.get(6)[1])+","+
//                        (tipExpression.get(7)[0]-tipExpression.get(7)[1])+","+
//                        (tipExpression.get(8)[0]-tipExpression.get(8)[1])+","+
//                        (tipExpression.get(9)[0]-tipExpression.get(9)[1])+","+
//                        tipNeo.get(0).intValue()+","+tipNeo.get(1).intValue()+","+tipNeo.get(2).intValue()+","+
//                        tipNeo.get(3).intValue()+","+tipNeo.get(4).intValue()+","+tipNeo.get(5).intValue()+","+
//                        tipNeo.get(6).intValue()+","+tipNeo.get(7).intValue()+","+tipNeo.get(8).intValue()+","+
//                        tipNeo.get(9).intValue()+"\n");
//                
//                        System.out.println("rep number: "+ rep + "    lamTime: " + (System.currentTimeMillis()-start)); start = System.currentTimeMillis();
//                }
                
                
//                     if(tipNeo.get(0)==0 && tipNeo.get(1)==0 && tipNeo.get(2)==0 && tipNeo.get(3)==0 && 
//                             tipNeo.get(4)==0 && tipNeo.get(5)==0 && tipNeo.get(6)==0 && tipNeo.get(7)==0){
//                        rawDatWriter.append(diffusionRate + ","  + pstar + "," + theta + ","+
//                            (tipExpression.get(0)[0]-tipExpression.get(0)[1])+","+
//                            (tipExpression.get(1)[0]-tipExpression.get(1)[1])+","+
//                            (tipExpression.get(2)[0]-tipExpression.get(2)[1])+","+
//                            (tipExpression.get(3)[0]-tipExpression.get(3)[1])+","+
//                            (tipExpression.get(4)[0]-tipExpression.get(4)[1])+","+
//                            (tipExpression.get(5)[0]-tipExpression.get(5)[1])+","+
//                            (tipExpression.get(6)[0]-tipExpression.get(6)[1])+","+
//                            (tipExpression.get(7)[0]-tipExpression.get(7)[1])+"\n");
//                
//                        System.out.println("rep number: "+ rep + "    lamTime: " + (System.currentTimeMillis()-start)); start = System.currentTimeMillis();
//                    }

//##############################################################################################################################################################################
              
//                outerDif[rep]=tipExpressionReplicates.get(0)[0]-tipExpressionReplicates.get(0)[1];
//                System.out.println(rep+ "   " + driftObject.getMuts()[1].intValue());
            
                tipExpression.clear();
                tipNeo.clear();
                driftObject.resetMuts();
                driftObject.resetPopNeo();
                driftObject.resetPopPsuedo();
                
                
               
            }// End loop for experiment reps
                
//           rawDatWriter.append(u + "," + theta + "," + omega + ","  + variance(outerDif) + "," + N + "," + brLength.get(0)+"\n");
           

            }// End outer loop (used for making replicas for finding omega/root(N) saturation)
            
        rawDatWriter.close();
        }
        catch(IOException L){L.printStackTrace();}
    } 
    
    
//###########################   END SIMULATION   ##############################################################    
    
    
    
// Some useful functions maybe 
    
    private static double mean(double[] dat){
        double sum = 0;
        for (int i = 0; i < dat.length; i++) {
            sum += dat[i];
        }
        return sum / dat.length;
    }
    private static double variance(double[] dat){
        double meanv = ExpressionDrift2.mean(dat);
        double temp = 0;
        for(double a : dat){
            temp += (meanv-a)*(meanv-a);
        }
        return temp/dat.length;
    }
    private static double[] aoverTotalcalculator (ArrayList<Double[]> x){
        double[] output = new double[x.size()];
        for(int i=0; i<x.size(); i++){
             output[i] = x.get(i)[0]/(x.get(i)[0]+x.get(i)[1]);
             if(x.get(i)[0]==0 && x.get(i)[1]==0){
                 output[i] = 0;
             }
        }
        return output;
    }
 
}
