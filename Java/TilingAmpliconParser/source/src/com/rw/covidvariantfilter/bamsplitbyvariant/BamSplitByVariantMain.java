/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.bamsplitbyvariant;

import com.rw.covidvariantfilter.Parameters;
import com.rw.covidvariantfilter.variants.OneVariantData;
import com.rw.globals.TerminalColors;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author raine
 */
public class BamSplitByVariantMain {

    public static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(BamSplitByVariantMain.class);

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
        Parameters.writeBams = cmd.hasOption("w");
        Parameters.variantTSV = cmd.getOptionValue("v");
        Parameters.variants = OneVariantData.parse(Parameters.variantTSV);
        if (cmd.hasOption("f")) {
            Parameters.minFractionOfMutationRequired = Optional.of(Float.valueOf(cmd.getOptionValue("f")));
        }
        Parameters.minMutationsInFragment = cmd.hasOption("m") ? Optional.of(Integer.valueOf(cmd.getOptionValue("m"))) : Optional.empty();
        Parameters.outputDir = cmd.getOptionValue("o");
        if (cmd.hasOption('i')) {
            String inBam = cmd.getOptionValue("i");
            if ((new File(inBam)).isFile()) {
                new BamSplitByVariant(inBam).parseBam();
            } else {
                try (Stream<Path> files = Files.walk(Path.of(inBam))) {
                    //key is bam file name, value map ; key variant, value parse stats
                    Map<String, Map<String, BamSplitByVariant.OneVariantParseStats>> stats = files.filter(Files::isRegularFile)
                            .parallel()
                            .filter(path -> path.toString().endsWith(".bam"))
                            .collect(Collectors.toConcurrentMap(
                                    k -> k.getFileName().toString(),
                                    v -> new BamSplitByVariant(v.toString()).parseBam()));

                    List<String> variants = stats.values().stream().flatMap(v -> v.keySet().stream()).distinct().sorted().toList();
                    List<String> samples = stats.keySet().stream().sorted().toList();
                    PrintWriter writer = new PrintWriter(new File(Parameters.outputDir, "Stats.tsv"));
                    writer.println(samples.stream().collect(Collectors.joining("\t\t", "\t", "")));
                    writer.print(IntStream.range(0, samples.size()).mapToObj(f -> "\tparsed\tmatching").collect(Collectors.joining()));
                    for (String variant : variants) {
                        writer.println();
                        writer.print(variant);
                        for (String sample : samples) {
                            BamSplitByVariant.OneVariantParseStats s = stats.get(sample).get(variant);
                            if (s == null) {
                                writer.print("\tN.D.\tN.D.");
                            } else {
                                writer.print("\t" + String.valueOf(s.parsed().get()) + "\t" + String.valueOf(s.matched().get()));
                            }
                        }
                    }
                    writer.println();
                    writer.close();
                } catch (IOException ex) {
                    LOGGER.fatal(ex);
                }
            }
        }
    }

    /**
     *
     * @return
     */
    private static Options cli_otions() {
        Options options = new Options();
        options.addOption(Option.builder("v")
                .required(true)
                .longOpt("variantData")
                .desc("tsv file with mutations defining variants")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("i")
                .required(true)
                .longOpt("inbam")
                .desc("input bam file or directory (will take all bam files recursively)")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("o")
                .required(true)
                .longOpt("outDir")
                .desc("output directory")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("m")
                .required(false)
                .longOpt("minmuts")
                .desc("minimal number of variant specific mutations to retain sam record")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("f")
                .required(false)
                .longOpt("fractionmutsrequired")
                .desc("fraction of mutations required to retain fragment")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("w")
                .required(false)
                .longOpt("writeBams")
                .desc("if set writes bam for each variant")
                .numberOfArgs(0)
                .build());
        return options;
    }
}
