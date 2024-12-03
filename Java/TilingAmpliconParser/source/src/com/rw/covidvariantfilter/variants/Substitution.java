/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.variants;

import htsjdk.samtools.SAMRecord;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.commons.lang3.tuple.ImmutablePair;

/**
 *
 * @author raine
 */
public class Substitution extends OneAlterationBase {

    public static Pattern regexGroup = Pattern.compile("([A-Za-z])(\\d+)([A-Za-z])");
    public static String testRegexString = ".?[a-zA-Z]\\d+[a-zA-Z]$";
   
    /*
       Matcher matcher = pattern.matcher(input);

        if (matcher.find()) {
            char firstChar = matcher.group(1).charAt(0);
            int number = Integer.parseInt(matcher.group(2));
            char secondChar = matcher.group(3).charAt(0);
     */
    public final char orgNuc;
    /**
     * list of possible nucs for this position, 
     */
    public final char altNuc;

    /**
     *
     * @param position
     * @param orgNuc
     * @param altNuc
     * @param score
     * @param aminoacidChange
     */
    public Substitution(int position, char orgNuc, char altNuc, int score, String aminoacidChange) {
        super(position, score, aminoacidChange);
        //type = score < 0 ? Type.SUBSTITUTION_ABSENCE : Type.SUBSTITUTION;
        this.orgNuc = orgNuc;
        this.altNuc = altNuc;
    }
    
    /**
     * factory method
     * @param mutstr
     * @return 
     */
      public static Substitution generate(String mutstr) {
        String[] s = mutstr.split("/");
        Matcher matcher = regexGroup.matcher(s[0]);
        if (matcher.find()) {
            char firstChar = matcher.group(1).charAt(0);
            int number = Integer.parseInt(matcher.group(2));
            //char[] altNuc = new char[1];
            char altNuc = matcher.group(3).charAt(0);
            
            return new Substitution(number, firstChar, altNuc, 1,  OneAlterationBase.extractAAmutString(mutstr));
        } 
        return null;
    }

    public static boolean matchesPattern(String s) {
        return s.matches(testRegexString);
    }

    
    @Override
    public Optional<OneAlterationBase> check(SAMRecord sam) {
        //int score = 0;
        final String read = sam.getReadString();
        int readpos = sam.getReadPositionAtReferencePosition(this.position);
        if (readpos > 0 && readpos < read.length()) { // readpos > read.length() - 1 this happens, WEIRD
            char readbase = read.charAt(readpos - 1);
            if (readbase == altNuc) 
                return Optional.of(this);
        }
        return Optional.empty();
    }
    
      @Override
    public boolean equals(Object obj) {
        if(obj instanceof Substitution == false)
            return false;
        Substitution s = (Substitution) obj;
        return orgNuc==s.orgNuc && altNuc == s.altNuc && position.intValue() == s.position.intValue(); 
    }

    @Override
    public int hashCode() {
        return new HashCodeBuilder(17, 37).
                append(orgNuc).
                append(altNuc).
                append(position).
                toHashCode();
    }

}
