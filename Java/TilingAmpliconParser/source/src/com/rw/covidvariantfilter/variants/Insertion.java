/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.variants;

import htsjdk.samtools.SAMRecord;
import java.util.Optional;
import org.apache.commons.lang3.builder.HashCodeBuilder;

/**
 *
 * @author raine
 */
public class Insertion extends OneAlterationBase {

    public final char[] insertion;

    /**
     *
     * @param position position where seq is inserted
     * @param insertion
     * @param absent
     * @param aminoacidChange
     */
    public Insertion(int position, String insertion, int score, String aminoacidChange) {
        super(position, score, aminoacidChange);
        // type = score < 0? Type.INSERTION_ABSENCE: Type.INSERTION;
        this.insertion = insertion.toCharArray();
    }

    @Override
    public Optional<OneAlterationBase> check(SAMRecord sam) {
        boolean passed = false;
        final String read = sam.getReadString();
        int readpos = sam.getReadPositionAtReferencePosition(this.position);
        if (readpos > 0 && sam.getReferencePositionAtReadPosition(readpos + this.insertion.length + 1) == this.position + 1) { //insertion length correct
            passed = true;
            for (int i = 0; i < this.insertion.length; i++) {
                if (read.charAt(readpos + i) != this.insertion[i]) //charat 0 based --> readpos corresponds to readpos +1
                {
                    passed = false;
                }
            }
        }
        if (passed) {
            return Optional.of(this);
        } else {
            return Optional.empty();
        }
    }

    @Override
    public boolean equals(Object obj) {
        boolean ok = false;
        if (obj instanceof Insertion ins) {
            ok = position.intValue() == ins.position.intValue();
            for (int i = 0; ok == true && i < insertion.length; i++) {
                ok &= insertion[i] == ins.insertion[i];
            }
        }
        return ok;
    }

    @Override
    public int hashCode() {
        return new HashCodeBuilder(17, 37).
                append(insertion).
                append(position).
                toHashCode();
    }
}
