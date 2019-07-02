/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package wrightfisher;

import java.util.ArrayList;

/**
 *
 * @author AmmonThompson
 */
public class BinarySearcher {
    public int finalX;
    public organism binarySearch(int up, int low, double randomDouble, ArrayList<Double> prob, ArrayList<organism> pop){
        int x = (up-low)/2+low+(up-low)%2;
        if(randomDouble<prob.get(0)){
            finalX=0;
        }
        else if(randomDouble>=prob.get(x-1) && randomDouble<prob.get(x)){
            finalX=x;
        }
        else{    
            if(prob.get(x)<randomDouble){
                binarySearch(up,x,randomDouble,prob,pop);
            }
            else{
                binarySearch(x,low,randomDouble,prob,pop);
            }
        }
        return new organism(pop.get(finalX).getMutations(), pop.get(finalX).getids(), pop.get(finalX).getExpression());
    }
}
