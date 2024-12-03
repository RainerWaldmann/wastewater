/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 */
package com.rw.mpileupparser;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
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
public class Main {
    final static public Parameters params = new Parameters();
    /**
     * 
     * @param args 
     */
    public static void main(String[] args) {
        Options options = cli_otions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException ex) {
            // Logger.getLogger(IlluminaOxfordBCUmiMerger.class.getName()).log(Level.SEVERE, null, ex);
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Usage ", options);
            System.out.println("Command line parsing error:\n  ");
            System.exit(1);
        }
        params.mPileupFileName = cmd.getOptionValue("i");
        params.outFile = cmd.getOptionValue("o");
        try {
            new MpileupFileParser().parse();
        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
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
                .longOpt("infile")
                .desc("input mpileup file , - reads from stdin")
                .numberOfArgs(1)
                .build());
        options.addOption(Option.builder("o")
                .required(false)
                .longOpt("out")
                .desc("output File")
                .numberOfArgs(1)
                .build());
        return options;
    }
}
