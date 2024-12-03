/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.lang3.tuple.ImmutablePair;


/**
 *
 * @author raine
 */
public class OneVariantData extends ArrayList<OneAlterationBase> {

    private static final org.apache.logging.log4j.Logger LOGGER = org.apache.logging.log4j.LogManager.getLogger(OneVariantData.class.getName());
    //public int min_alt;
    //public float min_fraction;


    /**
     * parses tsv that is used for variant definition in python script for
     * variant definitions
     *
     * @param in
     * @return
     */
    public static Map<String, OneVariantData> parse(String in) {
        try {
            //just get a list of lines from TSV
            final List<String> tsvLines = Files.lines(Paths.get(in)).
                    filter(x -> x.startsWith("#")==false). // skip lines starting with #
                    collect(Collectors.toList());
            //StringTokenizer st = new StringTokenizer(l.get(0),"\t");
            //find line number where mutations start in TSV
            final int lineIndexMutations = IntStream.range(0, tsvLines.size()).filter(d -> tsvLines.get(d).startsWith("mutations")).findFirst().getAsInt();
            //get list with variants, take firts line of TSV split by tab and skip first element           
            List<String> variants = Arrays.stream(tsvLines.get(0).split("\t")).skip(1).collect(Collectors.toList());
            //get a list with lists of mutations in each line of TSV
            final List<List<String>> mutationLines = tsvLines.stream().skip(lineIndexMutations).//skip lines til mutations are reached
                    map(d -> Arrays.stream(d.split("\t")).skip(1).collect(Collectors.toList())).
                    //filter(a -> a.isEmpty() == false && a.get(0).startsWith("#") == false).
                    collect(Collectors.toList());
/*
same index in list of variants and in each sublist of List<List<String>> mutationLines correspond to same variant
--> put stuff together in a Map where key is the variant name and value is an object of type OneVariantData which is essentially an ArrayList<OneAlterationBase> 
*/
            return IntStream.range(0, variants.size()).//get stream of indices,
                    filter(d -> variants.get(d).startsWith("#") == false). /*exclude variants commented wit '#'*/
                    
                    mapToObj(h ->           //assemble data for one variant
                            new ImmutablePair<String,OneVariantData>(variants.get(h), // left of Pair is variant name 
                            mutationLines.stream().               /* right is OneVariantData, to get it stream all lines with mutations and get element for index h*/      
                            filter(e -> e.size() > h).     /*assure that there is something at this index*/
                            map(t -> t.get(h).replaceAll("\\s+", "")). // retrieve corresponding mutation and remove whitespace
                            filter(y -> y.length() > 0 && y.startsWith("#") == false). //filter empty fields and commented mutations
                            map(r ->{ // map mutation string to ? extends OneAlterationBase
                            OneAlterationBase oneAlt = null;
                            if (Substitution.matchesPattern(r.split("/")[0])) {
                                oneAlt = Substitution.generate(r);
                            } else if (Deletion.matchesPattern(r.split("/")[0])) {
                                oneAlt = Deletion.generate(r);
                            } else { //TODO treat insertion !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                LOGGER.error("Could not parse: " + r);
                            }
                            return oneAlt;  
                         }).
                        filter(f -> f!= null). //is stream of OneAlterationBase here
                        collect(Collectors.toCollection(OneVariantData::new))) // generate OneVariantData Object for this variant
                    )
                    
                    .filter(m -> m.left.length()>0 /*empty variant name*/&& m.left.startsWith("#") == false /*commented variant name*/&& m.right.size()>0 /*no mutations*/) // do some filtering 
                    .collect(Collectors.toMap(k -> k.left, v -> v.right)); // transform into map where key is left of pair and value is right of pair
            
        } catch (IOException ex) {
            LOGGER.fatal("PARSING OF VARIANT INFO TSV FAILED" + ex);
        }
        return null;
    }
    

}




//  SparkSession spark = SparkSession.builder()
//    .master("local[*]")
//    .appName("SparkTest")
//    .getOrCreate();
//      Dataset<Row> csvDataToDF = spark.read()
//     .option("header", true)
//     .option("delimiter", "\t")    
//    .csv(in);          
//
//        return null;
//
//      try {
//     Files.lines((new File(in)).toPath(), Charset.forName("UTF-8"))
//    .map(line -> line.split("\t"))
//    .map(f -> )
//    .collect(Collectors.toList());    
//     }catch (IOException e) {
//         
//     }

//            return variantMutations.entrySet().stream().map(dataForOneVariant -> new 
//                ImmutablePair<String,OneVariantData> (dataForOneVariant.getKey(),
//                       dataForOneVariant.getValue().stream().map(x -> {
//                            OneAlterationBase oneAlt = null;
//                            if (Substitution.matchesPattern(x.split("/")[0])) {
//                                oneAlt = Substitution.generate(x);
//                            } else if (Deletion.matchesPattern(x.split("/")[0])) {
//                                oneAlt = Deletion.generate(x);
//                            } else {
//                                LOGGER.error("Could not parse: " + x);
//                            }
//                            return oneAlt;
//                       }).filter(f -> f!= null).collect(Collectors.toCollection(OneVariantData::new)))).
//                    collect(Collectors.toMap(k -> k.left, v -> v.right));

 

//try (DirectoryStream<Path> stream = Files.newDirectoryStream(new File(in).toPath(), "*.{c,h,cpp,hpp,java}")) {
//       for (Path entry: stream) {
//
//       }
//   } catch (DirectoryIteratorException ex) {
//       // I/O error encounted during the iteration, the cause is an IOException
//       throw ex.getCause();
//   }

        
//        CsvReadOptions.Builder builder = 
//	CsvReadOptions.builder(in)
//		.separator('\t')										// table is tab-delimited
//		.header(false)											// no header
//		.dateFormat("yyyy.MM.dd");  				// the date format to use. 
//
//CsvReadOptions options = builder.build();
//
//Table t1 = Table.read().usingOptions(options);

