/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.primertrim;

import com.rw.covidvariantfilter.Main;
import com.rw.covidvariantfilter.Parameters;
import com.rw.covidvariantfilter.bed.BedReader;
import com.rw.covidvariantfilter.primertrim.Primertrimmer;
import com.rw.globals.TerminalColors;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import static java.util.stream.Collectors.joining;
import java.util.stream.Stream;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author raine
 */
public class PrimerTrimmerMain {

    public void doJob(String[] args) {
        Options options = cli_otions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException ex) {
            // Logger.getLogger(IlluminaOxfordBCUmiMerger.class.getName()).log(Level.SEVERE, null, ex);
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Usage ", options);
            System.out.println(TerminalColors.RED_BOLD + "Command line parsing error:\n  " + TerminalColors.RESET + ex);
            System.exit(1);
        }
        Parameters.AMPLICON_EXTREMITY_FUZZYNESS = Integer.parseInt(cmd.getOptionValue("f"));
        Parameters.writeSplitBams = cmd.hasOption("s");
        Parameters.writeNonMatching = cmd.hasOption("n");
        Parameters.isSpikeInBam = cmd.hasOption("z");
        String inBam = cmd.getOptionValue("i");
        if(cmd.hasOption("o"))
            Parameters.outBamFile = Optional.of(cmd.getOptionValue('o'));
        //String outBam = cmd.getOptionValue("o");

        Parameters.ampliconData = BedReader.readInFiles(cmd.getOptionValue("b"), null);
        new Primertrimmer(inBam).trimPrimers();
        printTrimStats(inBam);
    }

    private void printTrimStats(String inBam) {
        String outNameRoot = inBam.substring(0, inBam.lastIndexOf("."));
        String outName = Parameters.isSpikeInBam ?  (outNameRoot + ".TrimStats_SpikeIn.tsv") : (outNameRoot + ".TrimStats.tsv");
        //System.out.println("****************  Writing Trimstats: " + outName + " Root=" + outNameRoot + " BAM=" + inBam);

        Function<Integer, String> tabGenerator = (count) -> Stream.generate(() -> "abc").limit(6 * count).collect(joining());
        AtomicInteger count = new AtomicInteger(0);
         try (BufferedWriter writer = new BufferedWriter(new FileWriter(outName))) {
            List<List<String>> lst = Parameters.ampliconData.entrySet().stream().map(h
                    -> Arrays.stream(
                            (h.getKey() + "\t\t\t\t" + "\n"
                                    + "pool 1" + "\t" + h.getValue().stream().filter(x -> x.getPoolNumber() == 1).mapToInt(v -> v.getTotalCount()).sum() + "\t\t\t" + "\n"
                                    + "pool 2" + "\t" + h.getValue().stream().filter(x -> x.getPoolNumber() == 2).mapToInt(v -> v.getTotalCount()).sum() + "\t\t\t" + "\n"
                                    + "PoolNo\t" + "SetNo\t" + "Begin\t" + "End\t" + "Count\n"
                                    + h.getValue().stream().map(w -> w.getCountPrintString()).collect(Collectors.joining("\n"))).split("\n")).collect(Collectors.toList())).//lists of lines here
                    sorted((o1, o2) -> Integer.compare(o2.size(), o1.size())).collect(Collectors.toList());
            List<Iterator<String>> iterlist = lst.stream().map(x -> x.iterator()).collect(Collectors.toList());
            while (iterlist.stream().anyMatch(x -> x.hasNext())) {
                String s = "";
                for (Iterator<String> it : iterlist) {
                    if (it.hasNext()) {
                        s += it.next() + "\t\t";
                    }
                }
                writer.write(s + "\n");
            }
        } catch (IOException ex) {
            Logger.getLogger(PrimerTrimmerMain.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @return
     */
    private static Options cli_otions() {
        Options options = new Options();
        options.addOption(Option.builder("i")
                .required(true)
                .longOpt("inbam")
                .desc("input bam file, - for stdin")
                .numberOfArgs(1)
                .build());
              options.addOption(Option.builder("o")
                .required(false)
                .longOpt("outbam")
                .desc("output bam file, defaults to infile + trimmed.bam")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("f")
                .required(true)
                .longOpt("fuzzyness")
                .desc("fuzzyness for amplicon ends / bed file match")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("b")
                .required(true)
                .longOpt("bedfiles")
                .desc("directory with bed files")
                .numberOfArgs(1)
                .build());
           options.addOption(Option.builder("s")
                .required(false)
                .longOpt("writesplitbams")
                .desc("write a bam for each bed file")
                .numberOfArgs(0)
                .build());
           options.addOption(Option.builder("z")
                .required(false)
                .longOpt("isSpikeInBam")
                .desc("when set, bam is a spikeIn filtered bam. Used to change Outfile extensions")
                .numberOfArgs(0)
                .build());
                      options.addOption(Option.builder("n")
                .required(false)
                .longOpt("writeNonMatching")
                .desc("write a bam for non matching")
                .numberOfArgs(0)
                .build());
//        options.addOption(Option.builder("o")
//                .required(true)
//                .longOpt("outbam")
//                .desc("output bam file")
//                .numberOfArgs(1)
//                .build());
        return options;
    }
}
