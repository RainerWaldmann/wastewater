/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter;

import com.rw.covidvariantfilter.bed.BedRecord;
import com.rw.covidvariantfilter.variants.OneVariantData;
import java.util.List;
import java.util.Map;
import java.util.Optional;

/**
 *
 * @author raine
 */
public class Parameters {
    /**
     * number of extra nucleotides to allow for assigning an amplicon to a bed record
     */
    public static int AMPLICON_EXTREMITY_FUZZYNESS = 15;
   public static final String PRIMERSET_SAMTAG = "PS";
   public static final String PRIMERPOOL_SAMTAG = "PP";// pool 1 or 2 typically
   public static final String PRIMERAMPLICONNUMBER_SAMTAG = "PN";// 
   public static final String PRIMERAMPLICONRANGE_SAMTAG = "PR";
   public static Optional<String> outBamFile = Optional.empty();
   public static Map<String, OneVariantData> variants;
   public static Map<String, List<BedRecord>> ampliconData;
   public  static Optional<Integer> minMutationsInFragment = Optional.empty();
   public static Optional<Float> minFractionOfMutationRequired = Optional.empty();
   //public boolean retainallminmuts =false;
   public static String variantTSV;
   public static boolean writeSplitBams;
   public static boolean writeNonMatching;
   public static boolean writeBams;
   public static String outputDir;
   /**
    * indicates whether bam is spike-in reads enriched -> changes file extensions
    */
   public static boolean isSpikeInBam = false;
}
