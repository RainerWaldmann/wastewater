/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter;

import com.rw.covidvariantfilter.primertrim.PrimerTrimmerMain;
import com.rw.covidvariantfilter.bamsplitbyvariant.BamSplitByVariantMain;
import java.util.Arrays;
//-v D:\\Documents\\NetBeansProjects\\CovidVariantFilter\\Pango  -r  "C:\Users\raine\OneDrive\Bureau\Wastewater\sarscov2.fa"  -i "C:\Users\raine\OneDrive\Bureau\Wastewater\Haliotis_22.12\barcode09.bam" -m 3 -f 0.8
/**
 *
 * @author raine
 */
public class Main {
    static final String TRIMPRIMERS = "trimprimers";
    static final String SPLITVARIANTS = "splitvariants";
    //static public final Parameters params = new Parameters();
    /**
     *
     * @param args run args
     */
    public static void main(String[] args) {
        try {
        if(args.length == 0)
            usage();
        if(args[0].equals(TRIMPRIMERS)){
            new PrimerTrimmerMain().doJob(Arrays.copyOfRange(args, 1, args.length, String[].class));
        } else if(args[0].equals(SPLITVARIANTS)){
            new BamSplitByVariantMain().doJob(Arrays.copyOfRange(args, 1, args.length, String[].class));
        } else
            usage();
        }catch (Exception e){
           e.printStackTrace();
           System.exit(1);
        }
    }
    
    private static void usage(){
        System.out.println("USAGE:");
        System.out.println(SPLITVARIANTS + "\t<options>" + "\t\tGenerate one Bam File per variant if variant defining mutations were found in an amplicon\n\t\t will only put variant specific amplicons into Bam");
        System.out.println(TRIMPRIMERS + "\t<options>" + "\t Trims primer specified in Bed files, optionally generates separate Bam file for each primer pool");
        System.exit(0);
    }
}
