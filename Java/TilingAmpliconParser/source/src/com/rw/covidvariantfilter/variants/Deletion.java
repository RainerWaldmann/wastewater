/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.variants;

import static com.rw.covidvariantfilter.variants.Substitution.testRegexString;
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
public class Deletion extends OneAlterationBase {

    public static Pattern regexGroup = Pattern.compile("(\\d+)del(\\d+)");
    public static String testRegexString = "^\\d+del\\d+$";

    public static boolean matchesPattern(String s) {
        return s.matches(testRegexString);
    }

    public static Deletion generate(String mutstr) {
        String[] s = mutstr.split("/");
        Matcher matcher = regexGroup.matcher(s[0]);
        if (matcher.find()) {
            int pos = Integer.parseInt(matcher.group(1));
            int len = Integer.parseInt(matcher.group(2));
            //String aaChange = s.length > 1? s[1]:"";        
            return new Deletion(pos, len, 1, OneAlterationBase.extractAAmutString(mutstr)); // TODO implement aaChange
        }
        return null;
    }

    final public int length;

    /**
     *
     * @param position is first position deleted
     * @param length
     * @param score
     * @param aminoacidChange
     */
    public Deletion(int position, int length, int score, String aminoacidChange) {
        super(position, score, aminoacidChange);
        // type = score < 0 ? Type.DELETION_ABSENCE:Type.DELETION;
        this.length = length;
    }

    @Override
    public Optional<OneAlterationBase> check(SAMRecord sam) {
        //start at -2 and check for deletion til pos+2 to account for mapping issues
        //TODO CHECK WHETHER OK SINCE DELETION POS IN JSON IS FIRST NUC DELETED
        int curretRefPosition = this.position - 2; //start at -2
        int readpos = sam.getReadPositionAtReferencePosition(curretRefPosition);
        if (readpos > 0) {
            for (int pos = ++readpos; pos < readpos + 5; pos++) {
                int newRefPos = sam.getReferencePositionAtReadPosition(pos);
                if (newRefPos != 0) {
                    if (newRefPos == curretRefPosition + 1 + this.length) {
                        return Optional.of(this);
                    }
                    curretRefPosition = newRefPos;
                }
            }
        }
        return Optional.empty();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Deletion == false) {
            return false;
        }
        Deletion s = (Deletion) obj;
        return length == s.length && position.intValue() == s.position.intValue();
    }

    @Override
    public int hashCode() {
        return new HashCodeBuilder(17, 37).
                append(length).
                append(position).
                toHashCode();
    }
}
