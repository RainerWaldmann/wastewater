/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.mpileupparser.modifications;



/**
 *
 * @author raine
 */
 public class Insertion extends ModificationsBase {

        final String insertion;

        public Insertion(int position, String insertion, boolean antisense, char qv) {
            super(position, "+" + insertion, antisense, qv);
            this.insertion = insertion;
        }

        @Override
        public int hashCode() {
            return insertion.hashCode() ^ super.hashCode();
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof Insertion == false) {
                return false;
            }
            return ((Insertion) obj).insertion.equals(insertion) && super.equals(obj);
        }
    }