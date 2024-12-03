/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.mpileupparser.modifications;


/**
 *
 * @author raine
 */
    public class Deletion extends ModificationsBase {

        final int length;

        public Deletion(int position, int length, boolean antisense, char qv) {
            super(position, "-" + "N".repeat(length), antisense, qv);
            this.length = length;
        }

        @Override
        public int hashCode() {
            return Integer.valueOf(length).hashCode() ^ super.hashCode();
        }

        @Override
        public boolean equals(Object obj) {
            if (obj instanceof Deletion == false) {
                return false;
            }
            return ((Deletion) obj).length == length && super.equals(obj);
        }
    }
