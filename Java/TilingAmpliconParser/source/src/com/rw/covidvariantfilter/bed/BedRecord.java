package com.rw.covidvariantfilter.bed;

import com.github.lindenb.jvarkit.tools.pcr.ReadClipper;
import com.rw.covidvariantfilter.Parameters;
import com.rw.covidvariantfilter.variants.Deletion;
import com.rw.covidvariantfilter.variants.Insertion;
import com.rw.covidvariantfilter.variants.OneAlterationBase;
import com.rw.covidvariantfilter.variants.OneVariantData;
import com.rw.covidvariantfilter.variants.Substitution;
import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.SimpleFeature;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import org.apache.commons.lang3.tuple.ImmutablePair;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author raine
 */
public class BedRecord { // TODO change to implement Locatable
    /**
     * used to speed up pow2 
     */
    static final int[] PRIME_NUMBERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 
        67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,127,131,137,139,149,151,157};
//    static {
//      POWERS_OF_2 = new int[31];
//      POWERS_OF_2[0] = 0;
//      for(int i = 1; i <31;i++)
//         POWERS_OF_2[i] =POWERS_OF_2[i-1]*2; 
//    }
    private final Integer fwdPrimerStart;
    private final Integer fwdPrimerEnd;
    private final Integer revPrimerStart;
    private final Integer revPrimerEnd;
    private final Integer setNumber;
    /**
     * same as contig in htsjdk locatable
     */
    private final String contig;
    /**
     * primer pool (typically one or two)
     */
    private final int poolNumber;

    public int getPoolNumber() {
        return poolNumber;
    }
    
    /**
     * key variant name, value list of alterations
     */
    private final Optional<Map<String, List<OneAlterationBase>>> variantsData;
    /**
     * key holds variant, value holds map of key:mutFoundPattern, value: Nrecords
     * found for this mutcount
     */
    private Optional<Map<String, Map<MutFoundPattern, AtomicInteger>>> variantCounts = Optional.empty();
    /**
     * total SAM records corresponding to PCR fragment
     */
    private final AtomicInteger totalCount = new AtomicInteger(0);

    public Integer getTotalCount() {
        return totalCount.get();
    }
    
    //private final VariantCount variantCount;
/**
 * 
 * @param fwdPrimerStart
 * @param fwdPrimerEnd
 * @param revPrimerStart
 * @param revPrimerEnd
 * @param setNumber
     * @param poolNumber
 * @param variants
 * @param contig seee htsjdk Locatable contig
 */
    public BedRecord(Integer fwdPrimerStart, Integer fwdPrimerEnd, Integer revPrimerStart, Integer revPrimerEnd, Integer setNumber, int poolNumber,
            final Map<String, OneVariantData> variants, final String contig) {
        this.fwdPrimerStart = fwdPrimerStart;
        this.fwdPrimerEnd = fwdPrimerEnd;
        this.revPrimerStart = revPrimerStart;
        this.revPrimerEnd = revPrimerEnd;
        this.setNumber = setNumber;
        this.poolNumber = poolNumber;
        this.variantsData = getVariantDataForAmplicon(variants);
        this.contig = contig;
    }

    /**
     * Gets variants that have mutations on this amplicon
     *
     * @param variants
     * @return map key is variant , value is list of alterations on this PCR
     * fragment
     */
    private Optional<Map<String, List<OneAlterationBase>>> getVariantDataForAmplicon(final Map<String, OneVariantData> variants) {
        if(variants == null)
            return Optional.empty();
        Map<String, List<OneAlterationBase>> m = variants.entrySet().stream().
                map((e) -> {
                    List<OneAlterationBase> lst = e.getValue().stream().filter((x) -> x.position >= this.fwdPrimerStart && x.position <= this.revPrimerEnd).
                            map((v) ->{
                                v.isOnPrimer = v.position <= this.revPrimerEnd || v.position >= this.revPrimerStart;
                                return v;
                            }).        
                            collect(Collectors.toList());
                    return new ImmutablePair<String, List<OneAlterationBase>>(e.getKey(), lst);
                }).filter(f -> f.right.isEmpty() == false).collect(Collectors.toMap(ImmutablePair::getLeft, ImmutablePair::getRight));
        return m.isEmpty() ? Optional.empty() : Optional.of(m);
    }
    
    /**
     * 
     * @param sam
     * @return pair with distance at 5p and 3p
     */
    public ImmutablePair<Integer,Integer> getDistanceFromAmplicon(SAMRecord sam) {
      return new ImmutablePair<>(sam.getAlignmentStart() - fwdPrimerStart, sam.getAlignmentEnd() - revPrimerEnd);
    }
    /**
     * checks whether SAMrecord matches amplicon
     *
     * @param sam
     * @return
     */
    private Boolean matchesAmplicon(SAMRecord sam) {
        final int refBegin = sam.getAlignmentStart();
        final int refEnd = sam.getAlignmentEnd();
        return refBegin >= fwdPrimerStart - Parameters.AMPLICON_EXTREMITY_FUZZYNESS && refBegin <= fwdPrimerStart + Parameters.AMPLICON_EXTREMITY_FUZZYNESS
                && refEnd >= revPrimerEnd - Parameters.AMPLICON_EXTREMITY_FUZZYNESS && refEnd <= revPrimerEnd + Parameters.AMPLICON_EXTREMITY_FUZZYNESS;
    }
    
    /**
     * checks whether SAMrecord matches amplicon
     *
     * @param sam
     * @param countIt indicates whether should be counted to avoid double counting 
     * @param primerSet 
     */
    public void primerTrimAmplicon(SAMRecord sam, boolean countIt, final String primerSet) {
            if(countIt)
            totalCount.incrementAndGet(); //TODO potential double addition between addSamRecord and here
          ReadClipper.clip(sam, new SimpleFeature(this.contig, this.fwdPrimerEnd+1, this.revPrimerStart-1));
          //sam.setAttribute(Parameters.PRIMERSET_SAMTAG, primerSet + "_p" + this.poolNumber + "_" + this.setNumber + "_" + this.fwdPrimerStart + "-" + this.revPrimerEnd);
          sam.setAttribute(Parameters.PRIMERSET_SAMTAG, primerSet);
          sam.setAttribute(Parameters.PRIMERPOOL_SAMTAG, this.poolNumber);
          sam.setAttribute(Parameters.PRIMERAMPLICONNUMBER_SAMTAG, this.setNumber);
          sam.setAttribute(Parameters.PRIMERAMPLICONRANGE_SAMTAG, this.fwdPrimerStart + "-" + this.revPrimerEnd);
    }
    


//    /**
//     * Parses and adds Sam record
//     *
//     * @param sam
//     * @param trimIt
//     * @param primerSet
//     */
//    public void addSamRecord(SAMRecord sam, boolean trimIt, final String primerSet) {
//        boolean matches = trimIt? primerTrimAmplicon(sam, false,primerSet) : matchesAmplicon(sam);
//        if (matches == false) {
//            return;
//        }
//        totalCount.incrementAndGet();//TODO potential double addition between primerTrimAmplicon and here
//        if (this.variantsData.isEmpty()) {
//            return;
//        }
//        String read = sam.getReadString();
//        for (Map.Entry<String, List<OneAlterationBase>> onevariant : this.variantsData.get().entrySet()) {
//            if((Main.params.minMutationsInFragment.isPresent() && onevariant.getValue().size() < Main.params.minMutationsInFragment.get()))
//                continue;
//            MutFoundPattern mutPattern = new MutFoundPattern(onevariant.getValue().size());
//            for(int i = 0; i< onevariant.getValue().size();i++){
//                OneAlterationBase a = onevariant.getValue().get(i);               
//               mutPattern.pattern[i] = switch (a.type) {
//                    case SUBSTITUTION ->
//                        checkSubstitution(sam, (Substitution) a, read);
//                    case DELETION ->
//                        checkDeletion(sam, (Deletion) a);
//                    case INSERTION ->
//                        checkInsertion(sam, (Insertion) a, read);
//               };
//            }
//            if(mutPattern.getMutCount() != 0){
//               if(variantCounts.isEmpty())
//                   variantCounts = Optional.of(new HashMap<>());
//               if(variantCounts.get().containsKey(onevariant.getKey()) == false)
//                   variantCounts.get().put(onevariant.getKey(),new HashMap<>());
//               if(variantCounts.get().get(onevariant.getKey()).containsKey(mutPattern)== false)
//                   variantCounts.get().get(onevariant.getKey()).put(mutPattern,new AtomicInteger(0));
//               variantCounts.get().get(onevariant.getKey()).get(mutPattern).incrementAndGet();
//            }
//        }
//                    ///!!!!!!!!!!!!!!!!!!!!!!!  TODO  CODE !!!!!!!!!!!!!!!!!!!!!
//    }

    public Integer getFwdPrimerStart() {
        return fwdPrimerStart;
    }

    public Integer getFwdPrimerEnd() {
        return fwdPrimerEnd;
    }

    public Integer getRevPrimerStart() {
        return revPrimerStart;
    }

    public Integer getRevPrimerEnd() {
        return revPrimerEnd;
    }

    public Integer getSetNumber() {
        return setNumber;
    }

    public Integer getAmpliconStart() {
        return fwdPrimerStart;
    }

    public Integer getAmpliconEnd() {
        return revPrimerEnd;
    }

     /**
     *
     * @param sam
     * @param substitution
     * @param read
     * @return score, negative if change should not be in variant, 0 if not found
     */
    private short checkSubstitution(SAMRecord sam, Substitution substitution, String read) {
        short retval = 0;
        int readpos = sam.getReadPositionAtReferencePosition(substitution.position);
        if (readpos < 1 || readpos > read.length() - 1) { // readpos > read.length() - 1 this happens, WEIRD
            retval = 0;
        } else {
            char readbase = read.charAt(readpos - 1);
            if(substitution.altNuc == readbase)
               retval = substitution.score.shortValue() ;  
//            for (char c : substitution.altNuc) {
//                if (c == readbase) {
//                    retval = substitution.score.shortValue() ;
//                    break;
//                }
//            }
        }
        return retval;
    }

    /**
     *
     * @param sam
     * @param deletion
     * @return score, negative if change should not be in variant, 0 if not found
     */
    private short checkDeletion(SAMRecord sam, Deletion deletion) {
        //start at -2 and check for deletion til pos+2 to account for mapping issues
        //TODO CHECK WHETHER OK SINCE DELETION POS IN JSON IS FIRST NUC DELETED
       short retval = 0;
        int curretRefPosition = deletion.position - 2; //start at -2
        int readpos = sam.getReadPositionAtReferencePosition(curretRefPosition);
        if (readpos < 1) {
            retval = 0;
        } else {
            for (int pos = ++readpos; pos < readpos + 5; pos++) {
                int newRefPos = sam.getReferencePositionAtReadPosition(pos);
                if (newRefPos != 0) {
                    if (newRefPos == curretRefPosition + 1 + deletion.length) {
                        retval = deletion.score.shortValue();
                        break;
                    }
                    curretRefPosition = newRefPos;
                }
            }
        }
        return retval;
    }

    /**
     *
     * @param sam
     * @param insertion
     * @param read
     * @return score, negative if change should not be in variant, 0 if not found
     */
    private short checkInsertion(SAMRecord sam, Insertion insertion, String read) {
        short retval = 0;
        int readpos = sam.getReadPositionAtReferencePosition(insertion.position);
        if (readpos < 1) {
            retval = 0;
        } else if (sam.getReferencePositionAtReadPosition(readpos + insertion.insertion.length + 1) == insertion.position + 1) { //insertion length correct
            retval = insertion.score.shortValue();
            for (int i = 0; i < insertion.insertion.length; i++) {
                if (read.charAt(readpos + i) != insertion.insertion[i]) //charat 0 based --> readpos corresponds to readpos +1
                {
                    retval = 0;
                }
            }
        }
        return retval;
    }
    /**
     * get count string
     * @return 
     */
    public String getCountPrintString(){
        return this.poolNumber +"\t" + this.setNumber   +"\t" + this.fwdPrimerStart + "\t" + revPrimerEnd + "\t" + totalCount.get();
    }
    
//    /**
//     * holds one alteration associated with one variant
//     */
//    public class OneAlteration extends OneAlterationBase {
//
//        /**
//         * indicates whether it is on primer
//         */
//        final boolean isOnPrimer;
//
//        /**
//         *
//         * @param b
//         * @param isOnPrimer
//         */
//        public OneAlteration(OneAlterationBase b, Boolean isOnPrimer) {
//            super(b);
//            this.isOnPrimer = false;
//        }
//    }
///////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * pattern of mutations found for this amplicon and a given variant
     * holds scores for each found position
     */
    private class MutFoundPattern{
        /**
         * scores for each alteration, 0 if not found
         */
      final Short[] pattern; 

        private MutFoundPattern(int maxMutCount) {
            this.pattern = new Short[maxMutCount];
            Arrays.fill(this.pattern, 0);
        }
/**
 * neg score adds -1, pos score adds 1
 * @return 
 */
        private int getMutCount(){
           return (int) Arrays.stream(pattern).mapToInt(c -> c>0? 1: c==0? 0:-1).sum();
        }
        
        @Override
        public boolean equals(Object obj) {
            if(obj instanceof MutFoundPattern == false)
                return false;
            MutFoundPattern m = (MutFoundPattern)obj;
            return Arrays.equals(this.pattern, m.pattern);
        }

        @Override
        public int hashCode() {
            int hash = 0;
            if(pattern.length <= PRIME_NUMBERS.length){
            for(int i = 0; i < pattern.length; i++)
              hash += pattern[i]* PRIME_NUMBERS[i]* i;  
            return hash;
            } else
                return Arrays.hashCode(pattern);
        }
        
    }
    /**
     *
     */
//    public class VariantCount{
//        /**
//         * n mutations expected for variant
//         */
//        public final Integer expectedMutations;
//        /**
//         * key is n mutations of expected found, value is n SAM records
//         */
//        public final Map<Integer,AtomicInteger> foundRecords;
//        public VariantCount(Integer expectedMutations) {
//            this.expectedMutations = expectedMutations;
//            this.foundRecords = new TreeMap<>();
//        }
//        
//        public void increment(int mutcount){
//            foundRecords.putIfAbsent(mutcount, new AtomicInteger(0)).incrementAndGet();           
//        }
//    }
}
