/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wrightfisher;

/**
 *
 * @author AmmonThompson
 */
import java.util.Random;
import java.lang.Math;
import java.util.ArrayList;
import java.lang.Double;
import java.lang.Integer;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

public class WrightFisher {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        File freqFile = new File("freq_dat.csv");
        try{
//        FileWriter freqWriter = new FileWriter(freqFile);
        ArrayList<organism> population = new ArrayList<organism>();
//        for(int g = 1; g<=(new Integer(args[2]).intValue()); g++){
        for(int g = 1; g<=1; g++){
//        int reps = new Integer(args[0]).intValue();
        int reps = 1000;
        int[] fixedreps = new int[reps];
        int fixed=0;
//        int generations=new Integer(args[3]).intValue();
        int generations=100000;
        int popsize=20;
//        double u=new Double(args[1]).doubleValue()*(double)g; 
        double u=0.02*(double)g; 
        double sample;
        double totalFitness;
        for(int q =0; q<reps; q++){
            int id = 1;
            for(int i=0; i<popsize; i++){
                population.add(new organism());
            }
            Random randomN = new Random();
            System.out.println("rep: "+q);
            for(int k=0; k<generations; k++){
        // Mutation
                for(organism i : population){
                    if(randomN.nextDouble()<u){
                        i.mutate(id);
                        id++;
                    }
                }
        // Sampling prob
                totalFitness=0;
                for(organism i : population){
                    totalFitness=totalFitness+i.getFitness();
                }
                ArrayList<Double> probArray = new ArrayList();
                probArray.add(new Double(population.get(0).getFitness()/totalFitness));
                for(int i=1; i<popsize; i++){
                    probArray.add(new Double(population.get(i).getFitness()/totalFitness)+probArray.get(i-1).doubleValue());
                }
//         Next Generation 
//########## Fast algorithm for randomly sampling from small population <= 1000 ####################
                ArrayList<organism> temp = new ArrayList();
                for(int j=0; j<popsize; j++){
                    sample=randomN.nextDouble();
                    if(sample<probArray.get(0)){
                            temp.add(new organism(population.get(0).getMutations(), population.get(0).getids(), population.get(0).getExpression()));
                    }
                    else{
                        for(int i=1; i<popsize; i++){
                            if(sample>=probArray.get(i-1) && sample<probArray.get(i)){
                                temp.add(new organism(population.get(i).getMutations(), population.get(i).getids(), population.get(i).getExpression()));
                                break;
                            }
                        }
                    }
                }
                population=temp;
                
//########## Fast algorithm for randomly sampling from larger pops > 1000 ###########                
//                ArrayList<organism> temp = new ArrayList();
//                for(int m =0; m<popsize; m++){
//                    temp.add(m,new BinarySearcher().binarySearch(popsize-1,0, randomN.nextDouble(), probArray, population));
//                }
//                population=temp;
//########################################################################################                

//######## create allele frequency distribution ############################################
//                if(q==0 && k%1000==0){
//                    // CODE HERE
//                    freqWriter.append("generation_"+k+",");
//                        for(Double b : alleleFreq(population)){
//                            System.out.println("frequency: "+b.doubleValue());
//                                freqWriter.append(b.doubleValue()+",");
//                        }
//                        freqWriter.append("\n");
//                    System.out.println("generation: "+k);
//                }
//###########################################################################################                
            fixedreps[q]=fixedreps[q]+getFixed(population);
            }
//###### keep track of mutations #############################################################            
//            fixed=0;
//            boolean isfixed=true;
//            if(!population.get(0).getids().isEmpty()){
//                for(int mutindx=0; mutindx < population.get(0).getids().size(); mutindx++){
//                    isfixed=true;
//                    for(int popindx=1; popindx < population.size(); popindx++){
//                        if(population.get(popindx).getids().size() <= mutindx){isfixed=false;break;}
//                        if(population.get(popindx).getids().get(mutindx).intValue() != population.get(0).getids().get(mutindx).intValue()){
//                            isfixed=false;break;
//                        }
//                    }
//                    if(isfixed){fixed++;}
//                }
//            }
//            fixedreps[q]=fixed;

//############## Check the above is working ########################################################################
//            System.out.println("######### fixed: "+fixed);
//            for(int i=0;i<population.size();i++){
//                if(population.get(i).getids().isEmpty()){
//                    System.out.println(population.get(i).getids().isEmpty());
//                }
//                else{
//                    System.out.println(population.get(i).getids().get(4));
//                }
//            }
//            System.out.println("###############  ");
//            System.out.println();

//##############################################################################################################            
            population.clear();
        }
//            File data = new File(args[4]+"data"+g+".csv");
            File data = new File("output8\\data"+(g+31)+".csv");
            FileWriter datawriter = new FileWriter(data);
            datawriter.append(u+","+generations+","+new organism().getSelection()+","+popsize+"\n");
            for(int i : fixedreps){
                datawriter.append(i+"\n");
            }
            datawriter.close();
//            freqWriter.close();
        }
        }catch(IOException L){L.printStackTrace();}
    
    }
    private static ArrayList<Double> alleleFreq(ArrayList<organism> pop){
        // sort largest to smallest by idnumber
        ArrayList<organism> soda = new ArrayList();
        for(organism i : pop){soda.add(i);}
        ArrayList<Double> freqDist = new ArrayList();
        int j=0;
        while(!soda.isEmpty()){
            //put non-stupid code here
            freqDist.add(1.0);
            if(soda.size()==1){
                soda.clear();
            }
            else{
                int i=1;
                while(i<soda.size()){
                    if(soda.get(0).getids().get(soda.get(0).getids().size()-1).doubleValue()==soda.get(i).getids().get(soda.get(i).getids().size()-1).doubleValue()){
                        freqDist.set(j,new Double(freqDist.get(j).doubleValue()+1.0));
                        soda.remove(i);
                        i--;
                    }
                    i++;
                }
                soda.remove(0);
            }
            j++;
        }
        return(freqDist);
    }
    private static int getFixed(ArrayList<organism> pop){
        int fixed=0;
        boolean isfixed=true;
        if(!pop.get(0).getids().isEmpty()){
            for(int mutindx=0; mutindx < pop.get(0).getids().size(); mutindx++){
                isfixed=true;
                for(int popindx=1; popindx < pop.size(); popindx++){
                    if(pop.get(popindx).getids().size() <= mutindx){isfixed=false;break;}
                    if(pop.get(popindx).getids().get(mutindx).intValue() != pop.get(0).getids().get(mutindx).intValue()){
                        isfixed=false;break;
                    }
                }
                if(isfixed){fixed++;}
            }
        }
        
//        if(fixed!=0){System.out.println("fixed: "+fixed);}
        System.out.println("######### fixed: "+fixed);
        for(int i=0;i<pop.size();i++){
            if(pop.get(i).getids().isEmpty()){
                System.out.println(pop.get(i).getids().isEmpty());
            }
            else{
                for(Integer j : pop.get(i).getids()){
                    System.out.print(j.intValue()+" ");
                }
                System.out.println();
            }
        }
        System.out.println("###############  ");
        System.out.println();
        
        
        for(int i=fixed-1;i>=0;i--){
            for(organism j : pop){
                j.getids().remove(i);
            }
        }

        
        return fixed;
    }
}
