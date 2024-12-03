/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.mpileupparser.modifications;

import java.util.Optional;
import java.util.function.Function;

/**
 *
 * @author raine
 */
 public abstract class ModificationsBase {

     public static Function<Character, Integer> getQV = (qvChar) -> (int) qvChar - 32;
        final int position;
        final String mutString;

        public String getMutString() {
            return mutString;
        }

        int fwdcount = 0;
        int revcount = 0;
        int fwdQVsum = 0;
        int revQVsum = 0;

    public int getFwdcount() {
        return fwdcount;
    }

    public int getRevcount() {
        return revcount;
    }

        
        
        public Optional<Float> getFWDqv() {
            return fwdcount == 0 ? Optional.empty() : Optional.of((float) fwdQVsum / fwdcount);
        }

        public Optional<Float> getREVqv() {
            return revcount == 0 ? Optional.empty() : Optional.of((float) revQVsum / revcount);
        }

        public float getMeanQV() {
            float sumQV = 0;
            int n = 0;
            Optional<Float> qV = getFWDqv();
            if (qV.isPresent()) {
                n++;
                sumQV += qV.get();
            }
            qV = getREVqv();
            if (qV.isPresent()) {
                n++;
                sumQV += qV.get();
            }
            return sumQV / n;
        }

        public ModificationsBase(int position, String mutString, boolean antisense, char qvChar) {
            this.position = position;
            this.mutString = mutString;
            incrementCounts(antisense, qvChar);
        }

        public void incrementCounts(boolean antisense, char qual) {
            if (antisense) {
                revcount++;
                revQVsum += getQV.apply(qual);
            } else {
                fwdcount++;
                fwdQVsum += getQV.apply(qual);
            }
        }

        public int getCounts() {
            return fwdcount + revcount;
        }

        @Override
        public boolean equals(Object obj) {
            ModificationsBase m = (ModificationsBase) obj;
            return m.position == position;
        }

        @Override
        public int hashCode() {
            return Integer.valueOf(position).hashCode();
        }
    }