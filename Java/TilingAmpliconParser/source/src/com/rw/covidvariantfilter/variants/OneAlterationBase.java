/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.variants;

import htsjdk.samtools.SAMRecord;
import java.util.Optional;

/**
 *
 * @author raine
 */
public abstract class OneAlterationBase implements Comparable<OneAlterationBase>{
   
     /**
         * will extract AA mut string from: * T17859C/orf1b:Y1464Y_BA.2.10.1,XBB
         * @param s input such as * T17859C/orf1b:Y1464Y_BA.2.10.1,XBB
         * @return  aa mut string or empty string if no aa mut
         */
        public static String extractAAmutString(String s){
            String[] sp = s.split("/");
            if(sp.length > 1)
                return sp[1].split("_")[0];
         return "";  
        }
    
//    public static enum Type {
//        SUBSTITUTION("S"),
//        INSERTION("I"),
//        DELETION("D");
//               
//        private final String tag;
//
//        public String getTag() {
//            return tag;
//        }
//
//        private Type(String s){
//                this.tag = s;
//        }
//    }
//    
//    public final Type type;
        
        
    /**
 * position 1-based
 * deletion: position 1-based after which deletion starts
 * insertion: position 1-based after which insertion starts
 */
    public final Integer position;
    /**
     * indicates how important presence of this mutation is
     */
    public final Integer score;
    public Optional<String> aminoAcidChange;
    /**
     * currently used for association of amplicon bedrecords and mutation to indicate whether it is on primer
     */
    public boolean isOnPrimer;
      
    /**
     * 
     * @param position
     * @param score
     * @param aminoAcidChange
     * @param type 
     */
     protected OneAlterationBase(Integer position, int score, String aminoAcidChan){
        this.position = position;
        this.score = score;
        this.aminoAcidChange = (aminoAcidChan == null) ? Optional.empty() : Optional.of(aminoAcidChan);

     }
  
 /**
     * @param sam
     * @returns this object if found else Optional.empty
     */
    abstract public Optional<OneAlterationBase> check(SAMRecord sam);
     
    @Override
    public int compareTo(OneAlterationBase o){
      return position.compareTo(o.position);
    }
    
    
}
