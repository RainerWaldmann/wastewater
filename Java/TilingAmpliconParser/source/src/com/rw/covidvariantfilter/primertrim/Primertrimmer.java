/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.primertrim;

import com.rw.covidvariantfilter.Parameters;
import com.rw.covidvariantfilter.bed.BedRecord;
import com.rw.covidvariantfilter.bamsplitbyvariant.BamReader;
import com.rw.covidvariantfilter.bamsplitbyvariant.BamWriter;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;

/**
 *
 * @author raine
 */
public class Primertrimmer {
    private final String inBam;
    /**
     * key is variant value: records passed
     */
    private int nRecords = 0;
    private int nRecordsTrimmed = 0;
    private int nRecordsMultifragment = 0;

    public Primertrimmer(String inBam) {
        this.inBam = inBam;
    }

    public void trimPrimers() {
        BamReader bamreader = new BamReader(inBam);
        bamreader.start();
        SAMFileWriterFactory factory = new SAMFileWriterFactory().
                setUseAsyncIo(true);
        String outfileName = Parameters.outBamFile.isPresent()? Parameters.outBamFile.get() : inBam.substring(0, inBam.lastIndexOf(".")) + "Trimmed.bam";
        SAMFileWriter writer = factory.makeBAMWriter(bamreader.getHeader(), false,
                new File(outfileName));
        SAMFileWriter nonMatchingWriter = null;
        if (Parameters.writeNonMatching) {
            SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
            String n = outfileName.substring(0, outfileName.lastIndexOf(".")) + "NonMatching.bam";
            nonMatchingWriter = factory.makeBAMWriter(bamreader.getHeader(), true,
                    new File(n));
        }
        BamWriter splitBamwriter = null; //only used for split writing
        if (Parameters.writeSplitBams) {
            splitBamwriter = new BamWriter(/*Parameters.ampliconData.keySet(),*/ outfileName.substring(0, outfileName.lastIndexOf(".")), bamreader.getHeader(),false, "");
            splitBamwriter.start();
        }
        SAMRecord sam = bamreader.get();
        while (sam != null) {
            nRecords++;

            Optional<List<ImmutablePair<String,Integer>>> r = findMatchingAmplicons(sam);
            if (r.isPresent()) {
                nRecordsTrimmed++;
                writer.addAlignment(sam);
                if (splitBamwriter != null) {
                    splitBamwriter.addSamToQueue(r.get().get(0).getLeft(), sam);
                }
                if (r.get().size() > 1 && r.get().get(1).right - r.get().get(0).right >60) { // only count if really different fragment TODO check whether logical
                    nRecordsMultifragment++;
                }
            } else if (Parameters.writeNonMatching) {
                nonMatchingWriter.addAlignment(sam);
            }
            sam = bamreader.get();
        };
        writer.close();
        if (Parameters.writeNonMatching) {
            nonMatchingWriter.close();
        }
        if (Parameters.writeSplitBams) {
            splitBamwriter.addSamToQueue(null, null);
            try {
                splitBamwriter.join();
            } catch (InterruptedException ex) {
                Logger.getLogger(Primertrimmer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        System.out.println(nRecords + " records parsed,  " + nRecordsTrimmed + " records primertrimmed,  " + nRecordsMultifragment + " records matching multiple fragments");
    }


   /**
    * 
    * @param sam
    * @return Pair contains name of bed and sum of distance from predicted start and stop
    */
    private Optional<List<ImmutablePair<String,Integer>>> findMatchingAmplicons(final SAMRecord sam) {
        List<Triple<Integer, String, BedRecord>> matchingAmplicons = Parameters.ampliconData.entrySet().stream().
                map(oneSet -> oneSet.getValue().stream().//stream data for one
                        map((oneBed) -> {
                            Optional<Triple<Integer, String, BedRecord>> retval;// triple with <distancesum,name, bedrecord>
                            ImmutablePair<Integer, Integer> dist = oneBed.getDistanceFromAmplicon(sam); // pair with distance at 5p and 3p
                            boolean pass = Math.abs(dist.getLeft()) <= Parameters.AMPLICON_EXTREMITY_FUZZYNESS
                                    && Math.abs(dist.getRight()) <= Parameters.AMPLICON_EXTREMITY_FUZZYNESS;
                            if (pass) {
                                retval = Optional.of(new ImmutableTriple(dist.getLeft() + dist.getRight(), oneSet.getKey(), oneBed));
                            } else {
                                retval = Optional.empty();
                            }
                            return retval;
                })
                ). //is stream of streams here
                flatMap(x -> x).
                filter(Optional::isPresent).map(d -> d.get()).
                sorted((o1, o2) -> o1.getLeft().compareTo(o2.getLeft())). // best match, shortest distance from an amplicon defined in a bed record first
                collect(Collectors.toList());
        if (matchingAmplicons.isEmpty() == false) {
            Triple<Integer, String, BedRecord> bestmatch = matchingAmplicons.get(0);
            bestmatch.getRight().primerTrimAmplicon(sam, true, bestmatch.getMiddle());
            return Optional.of(matchingAmplicons.stream().map(c -> new ImmutablePair<>(c.getMiddle(),c.getLeft())).collect(Collectors.toList()));
        } else {
            return Optional.empty();
        }
//                map(oneBed -> oneBed.primerTrimAmplicon(sam, true, oneSet.getKey()) == true ? oneSet.getKey() : "")).
//                flatMap(x -> x).
//                filter(f -> f.isEmpty() == false).
//                collect(Collectors.toList());
//        return (int) Parameters.ampliconData.values().stream().flatMap(x -> x.stream()).
//                map(c -> c.primerTrimAmplicon(sam,true)).filter(f -> f == true).count();
    }
}
