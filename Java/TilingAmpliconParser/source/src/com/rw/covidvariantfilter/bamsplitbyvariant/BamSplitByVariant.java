/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.bamsplitbyvariant;

import com.rw.covidvariantfilter.Main;
import com.rw.covidvariantfilter.Parameters;
import com.rw.globals.Globals;
import com.rw.covidvariantfilter.variants.OneAlterationBase;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author raine
 */
public class BamSplitByVariant {

    public static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(BamSplitByVariant.class);
    private final String inBam;
    

    /**
     * used for stats inner maps contains stats for individual amplicons
     */
    static public record OneVariantParseStats(AtomicInteger parsed, AtomicInteger matched, Map<String, OneVariantParseStats> ampliconstats) {

        public OneVariantParseStats() {
            this(new AtomicInteger(0), new AtomicInteger(0), new HashMap<>());
        }
    }
    ;
    /**
     * key is variant value: value is left: records that can contain required
     * mutations right: records that pass
     */
    private final Map<String, OneVariantParseStats> stats;

    /**
     *
     * @param inBam
     */
    public BamSplitByVariant(String inBam) {
        this.inBam = inBam;
        //populate stats with an entry for each variant
        stats = Parameters.variants.keySet().stream().collect(Collectors.toMap(k -> k, v -> new OneVariantParseStats()));
    }

    public Map<String,OneVariantParseStats> parseBam() {
        int n_samsParsed = 0;        
        final SamReaderFactory srf = SamReaderFactory.makeDefault();
        srf.validationStringency(ValidationStringency.SILENT);
        srf.setUseAsyncIo(true);
        SamReader sr = inBam.equals("-") ? srf.open(SamInputResource.of(System.in)) : srf.open(new File(inBam));

        BamWriter bamwriter = null;
        if (Parameters.writeBams) {
            bamwriter = new BamWriter(FilenameUtils.normalize(Parameters.outputDir) + "/", sr.getFileHeader(), true, inBam.equals("-") ? "" : (new File(inBam)).getName().split("\\.ba")[0] + "_");
            bamwriter.start();
        }
        
        for (SAMRecord sam : sr) {
            n_samsParsed++;
            //BamWriter.SAMwithID samWithId = new BamWriter.SAMwithID(sam);
            List<OneVariantParseResult> ret = parseOne(sam);
            if (bamwriter != null) {
                for (OneVariantParseResult r : ret) {
                    bamwriter.addSamToQueue(r.variant, sam);
                }
            }
        }
        if (bamwriter != null) {
            bamwriter.addSamToQueue(null, null); // poison element to terminate BamWriter Thread
            try {
                bamwriter.join();
            } catch (InterruptedException ex) {
                LOGGER.error(ex);
            }
        }
        System.out.println("\nSample; " + this.inBam);
        System.out.println(n_samsParsed + " records parsed, (Minmuts per sam: " + Parameters.minMutationsInFragment.get() + " Minfraction: " + Parameters.minFractionOfMutationRequired.get() + ")");

        try {
            String filenameRoot = (new File(inBam)).getName();
            filenameRoot = filenameRoot.substring(0, filenameRoot.lastIndexOf("."));
            PrintWriter writer = new PrintWriter( new File(Parameters.outputDir , filenameRoot + "_Min" + Parameters.minMutationsInFragment.get() + "_MinFrac" + Parameters.minFractionOfMutationRequired.get() + ".variants.tsv"));
            writer.println("Variant\tTested\tMatches");
            stats.entrySet().stream().forEach(stat -> {
                System.out.println(stat.getKey() + ":\ttested: " + stat.getValue().parsed.get() + "\tmatched: " + stat.getValue().matched.get());
                writer.println(stat.getKey() + "\t" + stat.getValue().parsed.get() + "\t" + stat.getValue().matched.get()
                        + stat.getValue().ampliconstats.entrySet().stream().map(e -> e.getKey() + ":" + e.getValue().parsed.get() + "/" + e.getValue().matched.get()).collect(Collectors.joining("\t", "\t", "")));
                /*           if (stat.getValue().getRight().get() == 0) {
                (new File(bamwriter.getOutFileName(stat.getKey()))).delete();
           } */
            });
            writer.close();
        } catch (IOException ex) {
            Logger.getLogger(BamSplitByVariant.class.getName()).log(Level.SEVERE, null, ex);
        }

        //System.out.println(nRecords + " records parsed,  " + nRecordsTrimmed + " records primertrimmed");
        return this.stats;
    }

    /**
     *
     * @param sam
     * @return list of passing variants with info on matches ....
     */
    private List<OneVariantParseResult> parseOne(SAMRecord sam) {
        final int refBegin = sam.getAlignmentStart();
        final int refEnd = sam.getAlignmentEnd();
        final String ampliconSetName = sam.getStringAttribute(com.rw.covidvariantfilter.Parameters.PRIMERSET_SAMTAG);
        final Integer ampliconPool = sam.getIntegerAttribute(com.rw.covidvariantfilter.Parameters.PRIMERPOOL_SAMTAG);
        final Integer ampliconNumber = sam.getIntegerAttribute(com.rw.covidvariantfilter.Parameters.PRIMERAMPLICONNUMBER_SAMTAG);
        final String ampliconRange = sam.getStringAttribute(com.rw.covidvariantfilter.Parameters.PRIMERAMPLICONRANGE_SAMTAG);
        final String ampliconIdentifier = ampliconSetName == null ? null : (ampliconSetName+ "_p" + ampliconPool + "_" + ampliconNumber + "_" + ampliconRange);
        // Map<Integer, List<OneSamParseResult>> groupedResult
        List<OneVariantParseResult> result = Parameters.variants.entrySet().stream().
                    map((oneVariant) -> {  // map to OneVariantParseResult
                        //get all mutations of one variant that can be on sam
                        List<OneAlterationBase> possibleMutationsOnSam = oneVariant.getValue().stream().
                            filter(c -> c.position >= refBegin + (Globals.PRIMERS_TRIMMED ? 0 : 25)
                            && c.position <= refEnd - (Globals.PRIMERS_TRIMMED ? 0 : 25)).
                            toList();

                        if (possibleMutationsOnSam.isEmpty() == false && (Parameters.minMutationsInFragment.isEmpty() || possibleMutationsOnSam.size() >= Parameters.minMutationsInFragment.get())) {
                            //increment the number of sam records that were analyzed for this variant
                            final OneVariantParseStats onevarstats = stats.get(oneVariant.getKey());
                            onevarstats.parsed.incrementAndGet();
                            //add stat for total parsed for amplicon
                            if (ampliconIdentifier != null) {
                                if (onevarstats.ampliconstats.containsKey(ampliconIdentifier)) {
                                    onevarstats.ampliconstats.get(ampliconIdentifier).parsed.incrementAndGet();
                                } else {
                                    onevarstats.ampliconstats.put(ampliconIdentifier, new OneVariantParseStats(new AtomicInteger(1), new AtomicInteger(0), null)); //one parsed , 0 matched so far
                                }
                            }

                            OneVariantParseResult oneVariantParseResult = possibleMutationsOnSam.stream().map((mut) -> mut.check(sam)).
                                    filter(m -> m.isPresent()).
                                    map(Optional::get).
                                    collect(Collector.of(() -> new OneVariantParseResult(oneVariant.getKey(), possibleMutationsOnSam), //supplier
                                        (accumulated, newitem) -> accumulated.mutationsFound.add(newitem),  //accumulator
                                        (part1, part2) -> part1.add(part2) // combiner
                                    ));
                            return oneVariantParseResult;
                        } else { // return OneSamParseResult with no matches
                            return new OneVariantParseResult(oneVariant.getKey(), possibleMutationsOnSam);//no matches 
                    }
                }) // end of map function
                .filter((x) ->
                         //no min muts specified --> passes if all possible muts found even if it is just one
                        (x.mutationsFound.size() == x.mutsExpected && Parameters.minMutationsInFragment.isEmpty())
                        || // at least min muts in fragment and no min fraction defined
                        (Parameters.minMutationsInFragment.isPresent() && Parameters.minFractionOfMutationRequired.isEmpty() && x.mutationsFound.size() >= Parameters.minMutationsInFragment.get())
                        || //both min muts and min fraction defined -> at least min muts and at least fraction of possible muts
                        (Parameters.minMutationsInFragment.isPresent() && Parameters.minFractionOfMutationRequired.isPresent()
                        && x.mutationsFound.size() >= Parameters.minMutationsInFragment.get() && x.mutationsFound.size() >= x.mutsExpected * Parameters.minFractionOfMutationRequired.get()))
                .map(a -> { //just add to stats
                            OneVariantParseStats onevarstats = stats.get(a.variant);
                            onevarstats.matched.incrementAndGet();
                            if (ampliconIdentifier != null) {
                                onevarstats.ampliconstats.get(ampliconIdentifier).matched.incrementAndGet();
                            }
                            return a;
                        })
                .toList();//collect(Collectors.groupingBy(OneSamParseResult::getScore)); //group by nmutsfound

        return result;

    }

    /**
     *
     */
    private class OneVariantParseResult implements Comparable<OneVariantParseResult> {

        private final String variant;
        /**
         * n mutations expected in this sam record for this variant
         */
        private final Integer mutsExpected;
        /**
         * each mutation can have a score this is the max score when all
         * mutations were found max possible score
         */
        //private final Integer maxScore;

        private final List<OneAlterationBase> mutationsFound;

        /**
         *
         * @param variant
         * @param possibleMuts
         */
        private OneVariantParseResult(String variant, List<OneAlterationBase> possibleMuts) {
            this.variant = variant;
            mutationsFound = new ArrayList<>();
            this.mutsExpected = possibleMuts.size();//(int) possibleMuts.stream().filter(e -> e.score > 0).count();
            //this.maxScore = (int) possibleMuts.stream().mapToInt(s -> s.score).sum();
        }

        /**
         * just add the found mutations
         *
         * @param v
         * @return
         */
        private OneVariantParseResult add(OneVariantParseResult v) {
            this.mutationsFound.addAll(v.mutationsFound);
            //this.mutsFound += v.mutsFound;
            //this.matchingAaMutsFound += v.matchingAaMutsFound;
            //this.score += v.score;
            //this.maxScore += v.maxScore;
            return this;
        }

        public String getVariant() {
            return variant;
        }

        public Integer getMutsExpected() {
            return mutsExpected;
        }

        /**
         *
         * @return
         */
        public Integer getScore() {
            return mutationsFound.stream().mapToInt(x -> x.score).sum();
        }

        /**
         * will compare number of muts found
         *
         * @param o
         * @return
         */
        @Override
        public int compareTo(OneVariantParseResult o) {
            return Integer.compare(this.mutationsFound.size(), o.mutationsFound.size());
        }

    }
}

// https://www.baeldung.com/java-19-record-patterns
//                     OneSamParseResult oneSamParseResult = mutations.stream().map((mut)
//                                -> switch (mut) {
//                            case Substitution s ->
//                                    s.check(sam);
//                            case Deletion d ->
//                                d.check(sam);
//                            case Insertion i ->
//                                i.check(sam);
//                            default -> null;
//                        }).
