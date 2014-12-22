/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.zooma;

/**
 * Provide a single entry point to C zooma code.
 *
 * To update the native library header file if you change the API here:
 * <pre>
 * ant compile-internal
 * cd /path/to/rtg-zooma/
 * make RTG_CLASSES=/path/to/rtg/build/src jni-h
 * </pre>
 *
 * To build the native library <code>libZooma.so</code>, see <a href="http://intranet.nz.realtimegenomics.com/wiki/pmwiki.php?n=RTG.ZoomaTodo">The Wiki Page</a>:
 *
 */
public class NativeZooma {

  private static final boolean IS_ENABLED;
  static {
    boolean loaded = false;
    try {
      System.loadLibrary("Zooma");
      loaded = true;
    } catch (UnsatisfiedLinkError e) {
      //Diagnostic.developerLog(e);
    }
    IS_ENABLED = loaded;
  }

  static boolean isEnabled() {
    return IS_ENABLED;
  }

  /**
   * Build a reference index
   * @param indexFile file name for output index file
   * @param templateFasta file name containing the reference FASTA
   * @param include string prefix for sequences to be included
   * @param exclude string prefix for sequences to be excluded
   * @return non-zero on failure
   */
  public native int buildIndex(String indexFile, String templateFasta, String include, String exclude);

  /**
   * Map paired-end reads
   * @param indexFile file name of already built index file
   * @param leftFasta file name of left arm FASTQ
   * @param rightFasta file name of left arm FASTQ
   * @param outPrefix name used as output prefix
   * @param readGroup read group specification
   * @param minFragSize minimum fragment size
   * @param maxFragSize maximum fragment size
   * @param threads number of threads to use
   * @param inputChunk number of reads per input chunk
   * @param outputChunk number of reads per output chunk
   * @param quickCache if true, build and use quick cache
   * @return non-zero on failure
   */
  public native int mapReads(String indexFile, String leftFasta, String rightFasta, String outPrefix, String readGroup, int minFragSize, int maxFragSize, int threads, int inputChunk, int outputChunk, boolean quickCache);

  /**
   * Map paired-end reads
   * @param indexFile file name of already built index file
   * @param leftFasta file name of left arm FASTQ
   * @param rightFasta file name of left arm FASTQ
   * @param outPrefix name used as output prefix
   * @param minFragSize minimum fragment size
   * @param maxFragSize maximum fragment size
   * @param threads number of threads to use
   * @param inputChunk number of reads per input chunk
   * @param quickCache if true, build and use quick cache
   * @return non-zero on failure
   */
  public native int mapfReads(String indexFile, String leftFasta, String rightFasta, String outPrefix, int minFragSize, int maxFragSize, int threads, int inputChunk, boolean quickCache);

}
