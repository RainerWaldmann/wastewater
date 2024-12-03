/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.covidvariantfilter.bamsplitbyvariant;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import java.io.IOException;
import java.util.Deque;
import java.util.concurrent.BlockingDeque;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author raine
 */
public class BamReader extends Thread{
    
    private final BlockingDeque<SAMRecord> queue;
    private final  SamReader sr;
    private boolean done = false;
    
    public BamReader(String in) {
        queue = new LinkedBlockingDeque<>(1024);
        final SamReaderFactory srf = SamReaderFactory.makeDefault();
     srf.validationStringency(ValidationStringency.SILENT);
    srf.setUseAsyncIo(true);
     sr =  in.equals("-")? srf.open(SamInputResource.of(System.in)): srf.open(new File(in));
    }
    
    public SAMFileHeader getHeader(){
        return sr.getFileHeader();
    }
    
    
    public void run(){       
     sr.forEach(record -> {
         try { 
             queue.put(record);
         } catch (InterruptedException ex) {
             Logger.getLogger(BamReader.class.getName()).log(Level.SEVERE, null, ex);
         }
     });
        try {
            sr.close();
        } catch (IOException ex) {
            Logger.getLogger(BamReader.class.getName()).log(Level.SEVERE, null, ex);
        }
     done = true;
     System.out.println("Bam Reader DONE");
    }
    
    public boolean isDone(){
        return done;
    }
 /**
  * 
  * @return 
  */   
    public Deque<SAMRecord> getDeque(){
        return queue;
    }
    /**
     * 
     * @return 
     */
    public synchronized SAMRecord get(){
        if(done && queue.isEmpty())
            return null;
        SAMRecord retval = null;
        try {
            retval = queue.poll(1, TimeUnit.SECONDS);
        } catch (InterruptedException ex) {
            Logger.getLogger(BamReader.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            
        }
        return retval;
    }
}
