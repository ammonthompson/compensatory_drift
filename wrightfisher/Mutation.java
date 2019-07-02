/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wrightfisher;

/**
 *
 * @author AmmonThompson
 */
public class Mutation {
    protected int genStart;
    protected int genFix;
    protected double mutationsize;
    public Mutation(int startGen, double size){
        genStart=startGen;
        mutationsize=size;
    }
    public void fixed(int fixtime){
        genFix=fixtime;
    }
    public int getStart(){return genStart;}
    public int getFix(){return genFix;}
    public double getMutsize(){return mutationsize;}
}
