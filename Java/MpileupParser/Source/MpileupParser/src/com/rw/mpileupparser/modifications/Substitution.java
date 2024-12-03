/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.mpileupparser.modifications;



/**
 *
 * @author raine
 */
public class Substitution extends ModificationsBase {

        final char refBase;
        final char newBase;

        public Substitution(int position, char refBase, char newBase, boolean antisense, char qv) {
            super(position, String.valueOf(newBase), antisense, qv);
            this.refBase = refBase;
            this.newBase = newBase;
        }

        @Override
        public int hashCode() {
            return refBase ^ newBase ^ super.hashCode();
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof Substitution == false) {
                return false;
            }
            Substitution s = (Substitution) obj;
            return s.refBase == refBase && s.newBase == newBase && super.equals(obj);
        }
    }
