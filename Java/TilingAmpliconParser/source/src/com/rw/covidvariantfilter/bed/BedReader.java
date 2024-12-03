/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.bed;

import com.rw.covidvariantfilter.variants.OneVariantData;
import com.rw.globals.TerminalColors;
import com.rw.globals.ThrowingFunction;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Triple;

/**
 * to read bed files with amplicon info
 *
 * @author raine
 */
public class BedReader {

    public static  HashMap<String, List<BedRecord>> readInFiles(String inDirectory, Map<String, OneVariantData> variants) {
        final BEDCodec BEDcodec = new BEDCodec();
        ThrowingFunction<String, Stream<Path>> getPaths = y -> Files.walk(Paths.get(y));
        ThrowingFunction<File, ImmutablePair<String, BufferedReader>> getReader
                = y -> new ImmutablePair<>(y.getName(), new BufferedReader(new FileReader(y)));
        
        HashMap<String, List<BedRecord>> bedFilesData = null;
        try{
        bedFilesData = getPaths.apply(inDirectory).
                filter(Files::isRegularFile).
                filter((p) -> p.toString().endsWith(".bed")).
                map(t -> t.toFile()).
                map(getReader).
                map((reader) -> {
                    return new ImmutablePair<String, List<BEDFeature>>(reader.getLeft(), reader.getRight().lines().
                            filter(StringUtils::isNotBlank).
                            filter((p) -> p.startsWith("browser") == false && p.startsWith("track") == false).//header
                            map(s -> s.replaceAll(" ", "")). //remove whitespace
                            //peek(h -> System.out.println(h)).
                            map(l -> BEDcodec.decode(l)).
                            collect(Collectors.toList()));
                }).
                collect(Collectors.toMap(
                        ImmutablePair::getLeft,
                        v -> transformToBedRecords.apply(v.getRight(), variants),
                        (o1, o2) -> o1,
                        HashMap::new));
        } catch (Exception e){
          System.err.println(" ++++++++++++++++++++++   Error parsing bed files with Amplicon Data - check files +++++++++++++++++++++++++++++++++++++"); 
          throw (e);
        }
        return bedFilesData;
    }
    /**
     * returns ordered list of amplicons from list of BEDFeatures
     */
    private static BiFunction<List<BEDFeature>, Map<String, OneVariantData>, List<BedRecord>> transformToBedRecords = (bedfeatures, variants) -> {
        //key is pairnumber // value left:pairnumber, middle isForward,right beadfeature
        Map<Integer, List<Triple<Integer, Boolean, BEDFeature>>> features = bedfeatures.stream().map((q) -> {
            String[] nameParts = q.getName().split("_");
            Triple<Integer, Boolean, BEDFeature> t = Triple.of(
                    Integer.parseInt(nameParts[nameParts.length - 2]),//pair number 
                    nameParts[nameParts.length - 1].startsWith("LEFT"), q);
            return t;
        }).collect(Collectors.groupingBy(Triple::getLeft));
        return features.entrySet().stream().
                map((onePair) -> {
                    Optional<Triple<Integer, Boolean, BEDFeature>> fwd = onePair.getValue().stream().filter((e) -> e.getMiddle()).findFirst();
                    Optional<Triple<Integer, Boolean, BEDFeature>> rev = onePair.getValue().stream().filter((e) -> e.getMiddle() == false).findFirst();
                    if (fwd.isEmpty() || rev.isEmpty()) {
                        String reason = fwd.isEmpty() && rev.isEmpty() ? "Both empty" : fwd.isEmpty() ? rev.get().getRight().getName() + " No matching forward primer"
                                : fwd.get().getRight().getName() + " No matching reverse primer";
                        TerminalColors.printErrorAndExit("EXITING: inconsistent BED !!!!!!!!!!!!!!!!!!!!!!!\n + reason"
                        );
                    }
                    return new BedRecord(fwd.get().getRight().getStart(), fwd.get().getRight().getEnd(),
                            rev.get().getRight().getStart(), rev.get().getRight().getEnd(), fwd.get().getLeft(), (int) fwd.get().getRight().getScore() /*pool*/, variants, fwd.get().getRight().getContig());
                }).sorted((o1, o2) -> o1.getAmpliconStart().compareTo(o2.getAmpliconStart())).collect(Collectors.toList());
    };

    
}
