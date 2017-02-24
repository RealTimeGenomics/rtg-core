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

import java.io.FileDescriptor;

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
//      e.printStackTrace(System.err);
//      Diagnostic.developerLog(e);
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
   * @param word word size for hashes
   * @param step step size between start of hashes
   * @return non-zero on failure
   */
  public native int buildIndex(String indexFile, String templateFasta, String include, String exclude, int word, int step);

  /**
   * Map paired-end reads
   * @param indexFilename filename for reference index
   * @param leftFasta file name of left arm FASTQ
   * @param rightFasta file name of left arm FASTQ
   * @param outPrefix name used as output prefix
   * @param numTheads number of threads
   * @param outputChunksize number of reads per chunk
   * @param eScore alignment score threshold
   * @param minFragSize minimum fragment size
   * @param maxFragSize maximum fragment size
   * @param step step size (<code>MAX_INT</code> to use indexes word size)
   * @param readGroup read group specification
   * @param buildQuick performance option
   * @param buildCache performance option
   * @param verbose extra debug output
   * @param threadMerge merge files per thread
   * @param fileMerge merge alls file into one
   * @param compressionLevel level of compression for output
   * @param refNamesFile name of reference sequences that should get own name (some output modes only)
   * @param stderr handle to where standard error should be sent (unused)
   * @param iCacheGB size of input cache in gigabytes
   * @param oCacheGB size of output cache in gigabytes
   * @return non-zero on failure
   */
  public native int mapReads(String indexFilename, String leftFasta, String rightFasta, String outPrefix, int numTheads,
                              int outputChunksize, int eScore, int minFragSize, int maxFragSize, int step,
                              String readGroup, boolean buildQuick, boolean buildCache, boolean verbose, boolean threadMerge,
                              boolean fileMerge, int compressionLevel, String refNamesFile, FileDescriptor stderr, int iCacheGB,
                              int oCacheGB);

}
