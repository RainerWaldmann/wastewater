/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.bamsplitbyvariant;

import com.rw.covidvariantfilter.Main;
import com.rw.covidvariantfilter.Parameters;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.tuple.ImmutablePair;

/**
 *
 * @author raine
 */
public class BamWriter extends Thread {
    /**
     * generates unique SAM record ID
     */
     private static AtomicLong SAM_ID = new AtomicLong(0);
    /**
     * holds samrecord and Id to assure ordered writing of records
     */
    public static record SAMwithID(SAMRecord sam, Long id) {
        public SAMwithID(SAMRecord sam) {
            this(sam, SAM_ID.incrementAndGet());
        }
    
    }

    public static final org.apache.logging.log4j.Logger LOGGER = org.apache.logging.log4j.LogManager.getLogger(BamWriter.class);
    private final SAMFileHeader header;
    private final Map<String, SAMFileWriter> samWriters = new HashMap();
    private final String prefix;
    /**
     * just used to generate outfile names
     */
    private final String fileNameRoot;
    /**
     * pair of name of variant or primer set and sam record
     */
    private final BlockingQueue<ImmutablePair<String, SAMRecord>> pendingBams = new LinkedBlockingDeque<>(1024);

    /**
     *
     * @param fileNameRoot to get the root of outfile name
     * @param header need header for new bam file
     * @param createIndex if true will create bai index file
     * @param prefix prefix to add to bam name
     */
    public BamWriter(String fileNameRoot, SAMFileHeader header, boolean createIndex, String prefix) {
        this.header = header;
        this.fileNameRoot = fileNameRoot;
        SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(createIndex);
        this.prefix = prefix;
    }

    @Override
    public void run() {
        SAMFileWriterFactory factory = new SAMFileWriterFactory();

        while (true) {
            if (pendingBams.isEmpty() == false) {
                ImmutablePair<String, SAMRecord> sam = pendingBams.remove();
                if (sam != null) {
                     if (sam.left == null && sam.right == null) {//null pair received -> terminate
                    break;
                }
                    if (samWriters.containsKey(sam.left) == false) {
                        samWriters.put(sam.left,
                                factory.makeBAMWriter(header, false, new File(getOutFileName(prefix + sam.left))));
                    }
                    samWriters.get(sam.left).addAlignment(sam.right);
                }
            } else {
                try {
                    Thread.sleep(10);
                } catch (InterruptedException ex) {
                    Logger.getLogger(BamWriter.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        for(SAMFileWriter w : samWriters.values())
            w.close();
    }

    /**
     *
     * @param suffix variant name or primer pool name
     * @return
     */
    private String getOutFileName(String suffix) {
        return fileNameRoot + suffix + (Parameters.isSpikeInBam ? "_SpikeIn" : "") + ".bam";
    }

    /**
     * parameters null null to signal end
     * will block if no space on queue
     * @param variant
     * @param sam send null null to terminate
     */
    public void addSamToQueue(String variant, SAMRecord sam) {
        try {
            pendingBams.put(new ImmutablePair<>(variant, sam));
        } catch (InterruptedException ex) {
            Logger.getLogger(BamWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
