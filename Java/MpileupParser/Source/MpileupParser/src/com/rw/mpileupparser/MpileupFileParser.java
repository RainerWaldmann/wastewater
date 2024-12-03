/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.rw.mpileupparser;

import com.rw.mpileupparser.modifications.Deletion;
import com.rw.mpileupparser.modifications.Insertion;
import com.rw.mpileupparser.modifications.ModificationsBase;
import com.rw.mpileupparser.modifications.Substitution;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.lang3.tuple.ImmutablePair;

/**
 *
 * @author raine
 */
public class MpileupFileParser {

//    enum Pattern {
//        START("^."),
//        END("$"),
//        INSERTION("+"),
//        IDENTICAL(",");
//
//        private final String pattern;
//
//        private Pattern(String s) {
//            this.pattern = s;
//        }
//    }
        /**
     * same pattern for insertion or deletion
     */
    final Pattern insertionPattern = Pattern.compile("^[0-9]+[ACGTacgt]+");
    //final Pattern senseInserion = Pattern.compile("^[AGCT]*");
    //final Pattern antisenseInserion = Pattern.compile("^[agct]*");
    /**
     * set to true when reader done
     */
    //AtomicBoolean readerDone = new AtomicBoolean(false);
    /**
     * holds some lines of the mpileup file
     */
    //private final BlockingQueue<String> linesQueue;

    /**
     *
     */
    public MpileupFileParser() {
        //linesQueue = new LinkedBlockingQueue<>(10);
    }

    public void parse() throws IOException {
        Stream<String> instream = Main.params.mPileupFileName.equals("-") ?
                new BufferedReader(new InputStreamReader(System.in)).lines() : 
                Files.lines(Path.of(Main.params.mPileupFileName));
        Map<Integer,MpileUpOnePositionInfo> mods = instream.parallel().map( (line) -> {
        //while (readerDone.get() == false && linesQueue.isEmpty() == false) {
            String[] ll = line.split("\t");
            int pos = Integer.parseInt(ll[1]);
            char refBase = ll[2].charAt(0);
            //int nReads = Integer.parseInt(ll[3]);
            final MpileUpOnePositionInfo mpileUpOnePositionInfo = new MpileUpOnePositionInfo(refBase);
            String patternString = ll[4];
            String qvString = ll[5];
            for (int i = 0, nth_read = 0; i < patternString.length(); i++) {
                char c = patternString.charAt(i);
                switch (c) {
                    case '#' -> { //deletion , nothing to do
                        nth_read++; 
                        mpileUpOnePositionInfo.incrementTotalReadsOnPos();
                    }
                    case '*' -> { //deletion , nothing to do
                        nth_read++;
                        mpileUpOnePositionInfo.incrementTotalReadsOnPos();
                    }
                    /*New read*/
                    case '^' -> {
                        //currentSeqIds.add(nthentry, currentSamRecordId++);
                        i++;//skip next char which is mapping quality
                    }
                    case ',' -> {
                        mpileUpOnePositionInfo.addRefMatch(true, qvString.charAt(nth_read));
                        mpileUpOnePositionInfo.incrementTotalReadsOnPos();
                        nth_read++;
                    }
                    case '.' -> {
                        mpileUpOnePositionInfo.addRefMatch(false, qvString.charAt(nth_read));
                        mpileUpOnePositionInfo.incrementTotalReadsOnPos();
                        nth_read++;
                    }
                    case '$' -> {
                        //toRemoveNextRound.add(nthentry - 1);
                        //nthentry++;
                    }
                    case '+', '-' -> {  //plus or minus refers to read corresponding to previous identity, substitution....
                       String s = patternString.substring(i+1);
                       Matcher matcher = insertionPattern.matcher(s);
                       matcher.find();
                       String indelString = matcher.group();
                       String alphabets = indelString.replaceAll("[^a-zA-Z]", "");
                       String digits = indelString.replaceAll("[^0-9]", "");
                       //String[] parts = indelString.split("(?<=\\D)(?=\\d)"); //splits into digits and non digits 
                       int indellength = Integer.valueOf(digits);
                       String indel = indelString.substring(digits.length(), digits.length() + indellength); // could be followed by substitution which is also an AGCT 
                       i += digits.length() /*+ indellength */ + alphabets.length();
                       //no increase of nthentry since info on match of current base follows on line
                       boolean reverse = indel.equals(indel.toUpperCase()) == false;
                       ModificationsBase mod = c == '+' ? new Insertion(pos, indel.toUpperCase(), reverse, qvString.charAt(nth_read-1)) : new Deletion(pos, indellength, reverse,qvString.charAt(nth_read-1) );
                       mpileUpOnePositionInfo.addModification(mod,reverse, qvString.charAt(nth_read-1)); //nthentry was incremented in previous loop -> -1
                    }
                    case 'a','g','c','t' -> { //antisense substitution
                       Substitution s = new Substitution(pos, refBase, Character.toUpperCase(c), true,qvString.charAt(nth_read));
                       mpileUpOnePositionInfo.addModification(s,true,qvString.charAt(nth_read));
                       mpileUpOnePositionInfo.incrementTotalReadsOnPos();
                       nth_read++;                                             
                    }
                     case 'A','G','T','C' -> { //sense substitution
                       Substitution s = new Substitution(pos, refBase, c, false,qvString.charAt(nth_read));
                       mpileUpOnePositionInfo.addModification(s,false,qvString.charAt(nth_read));
                       mpileUpOnePositionInfo.incrementTotalReadsOnPos();
                       nth_read++;
                    }


           
                }
            }
//               if (toRemoveNextRound.isEmpty() == false) {
//                for (int i = toRemoveNextRound.size() - 1; i >= 0; i--) {
//                    currentSeqIds.remove((int)(toRemoveNextRound.get(i)));
//                }
//                toRemoveNextRound.clear();
//            }
               //modifications.put(pos,mpileUpOnePositionInfo);
               return new ImmutablePair<Integer,MpileUpOnePositionInfo>(pos, mpileUpOnePositionInfo);
        }).collect(Collectors.toConcurrentMap(k -> k.getLeft(), v-> v.getRight()));
        Map<Integer,MpileUpOnePositionInfo> sortedMods = new TreeMap<Integer,MpileUpOnePositionInfo>(mods);
        writeTSV(sortedMods);
    }
    
    private void writeTSV(Map<Integer,MpileUpOnePositionInfo> sortedMods){
        DecimalFormat df = new DecimalFormat("#.0000");
        try(FileWriter writer = new FileWriter(Main.params.outFile)){
        writer.write("REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\n");       
       sortedMods.entrySet().stream().filter(x -> x.getValue().modifications != null).forEach((entry) -> {
           int pos = entry.getKey();
            MpileUpOnePositionInfo v = entry.getValue();
          entry.getValue().modifications.values().stream().filter(g -> (float)g.getCounts()/v.getTotalReadsOnPos() >= Parameters.minALTfreqToReport).forEach( oneModEntry -> {           
               try {
                   writer.write("\t" + pos + "\t"  + v.refBase  + "\t"  +  oneModEntry.getMutString() + "\t" + v.getRefDepth()  + "\t" +
                           v.getMatchesRefReverse()  + "\t" +  df.format(v.getRefQV()) + "\t" + oneModEntry.getCounts() + "\t" +  
                           oneModEntry.getRevcount() + "\t" +  oneModEntry.getMeanQV()  + "\t" +  df.format((float)oneModEntry.getCounts()/v.getTotalReadsOnPos()) + "\t" +  v.getTotalReadsOnPos() + "\n");
               } catch (IOException ex) {
                   Logger.getLogger(MpileupFileParser.class.getName()).log(Level.SEVERE, null, ex);
               }
          });
       });
       writer.close();
        } catch(IOException e){
            e.printStackTrace();
        }
    }

//    private class MpileupReader implements Runnable {
//
//        @Override
//        public void run() {
//            try (Stream<String> stream = Files.lines(Paths.get(Main.params.mPileupFileName))) {
//                stream.forEach(x -> {
//                    try {
//                        linesQueue.put(x);
//                    } catch (InterruptedException ex) {
//                        Logger.getLogger(MpileupFileParser.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                });
//            } catch (IOException e) {
//                System.err.format("Exception: %s%n", e);
//            }
//            readerDone.set(true);
//        }
//
//    }

}
