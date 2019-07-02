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
public class organism implements Cloneable{
    ArrayList<Double> mutations = new ArrayList();
    ArrayList<Integer> mutationid = new ArrayList();
    //ArrayList<Mutation> muts = new ArrayList();
    private double fitness;
    private double expression;
    private double omega=5;
    private Random randomN;
    public organism(ArrayList<Double> m, ArrayList<Integer> id, double e){
//        for(Double i : m){
//            mutations.add(new Double(i.doubleValue()));
//        }
        for(Integer i : id){
            mutationid.add(new Integer(i.intValue()));
        }
        expression=e;
        fitness = Math.exp(-0.5*Math.pow((expression-100)/omega,2));
//        fitness=0.01;
    }
    public organism(){
        fitness = 1;
        expression = 100;
    }
    public void mutate(int mutid){
        mutations.add(new Double(new Random().nextGaussian()));
        mutationid.add(mutid);
        expression=expression+mutations.get(mutations.size()-1).doubleValue();
        fitness = Math.exp(-0.5*Math.pow((expression-100)/omega,2));
//        fitness=0.01;
    }
    public double getFitness(){
        return(fitness);
    }
    public double getExpression(){
        return(expression);
    }
    public ArrayList<Double> getMutations(){
        return mutations;
    }
    public ArrayList<Integer> getids(){
        return mutationid;
    }
    public double getSelection(){
        return omega;
    }
    
}
