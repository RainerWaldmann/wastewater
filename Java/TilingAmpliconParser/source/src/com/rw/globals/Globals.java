/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.rw.globals;

import com.rw.globals.ThrowingConsumer;
import java.util.function.Consumer;

/**
 *
 * @author raine
 */
public class Globals {
   public final static String MUTATIONS_EXPECTEDINFRAGMENT_TAG = "NM";
   public final static String MUTATIONS_FOUNDINFRAGMENT_TAG = "FM";
   public final static String VARIANTS_MATCHING_TAG = "VA";
   public final static boolean PRIMERS_TRIMMED = true;
 
   
   /**
    * https://www.baeldung.com/java-lambda-exceptions
    * @param <T>
    * @param <E>
    * @param throwingConsumer
    * @param exceptionClass
    * @return 
    */
   public static <T, E extends Exception> Consumer<T> handlingConsumerWrapper(
  ThrowingConsumer<T, E> throwingConsumer, Class<E> exceptionClass) {
 
    return i -> {
        try {
            throwingConsumer.accept(i);
        } catch (Exception ex) {
            try {
                E exCast = exceptionClass.cast(ex);
                System.err.println(
                  "Exception occured : " + exCast.getMessage());
            } catch (ClassCastException ccEx) {
                throw new RuntimeException(ex);
            }
        }
    };
   }
}
