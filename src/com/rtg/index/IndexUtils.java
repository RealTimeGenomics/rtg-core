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
package com.rtg.index;

import com.rtg.index.params.CreateParams;
import com.rtg.index.params.ParamsUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Some utility functions used by index implementation.
 *
 */
public final class IndexUtils {

  private IndexUtils() { }

  /**
   * Create a new <code>Index</code> using the current default implementation.
   * @param indexParams holds all the values needed for constructing the index.
   * @param filter the index filter
   * @param threads ignored
   * @return the new <code>Index</code>
   */
  public static Index createIndex(final CreateParams indexParams, IndexFilterMethod filter, final int threads) {
    if (indexParams.compressHashes()) {
      Diagnostic.developerLog("Creating compressed index");
      return new IndexCompressed(indexParams, filter, threads);
    } else {
      Diagnostic.developerLog("Creating simple index");
      return new IndexSimple(indexParams, filter, threads);
    }
  }

  /**
   * Compute the total number of bytes of memory that the <code>IndexImplementation</code> will
   * require when it is constructed. Only includes terms that depend on size and ignores constant
   * terms.
   * @param createParams parameters that will be used to construct the <code>IndexImplementation</code>.
   * @return to the total bytes of memory.
   */
  public static long bytes(final CreateParams createParams) {
    long total = 0;
    total += createParams.hash().bytes();
    total += createParams.value().bytes();
    total += createParams.initialPosition().bytes();
    final HashBitHandle bitVector = createParams.bitVector();
    if (bitVector != null) {
      total += bitVector.bytes();
    }
    return total;
  }

  /**
   * Create a human readable description of the memory usage of an <code>IndexImplementation</code>
   * constructed using this.
   * @param createParams parameters that will be used to construct the <code>IndexImplementation</code>.
   * @return the memory description.
   */
  public static String memToString(final CreateParams createParams) {
    final StringBuilder sb = new StringBuilder();
    memToString(sb, createParams);
    return sb.toString();
  }

  /**
   * Create a human readable description of the memory usage of an <code>IndexImplementation</code>
   * constructed using this.
   * This version is intended for use in logging and conforms to a format that allows easy
   * extraction by text processing tools.
   * @param sb where to place the result.
   * @param createParams parameters that will be used to construct the <code>IndexImplementation</code>.
   */
  public static void memToString(final StringBuilder sb, final CreateParams createParams) {
    sb.append(ParamsUtils.memToString("Hash", createParams.hash().bytes()));
    sb.append(ParamsUtils.memToString("Value", createParams.value().bytes()));
    sb.append(ParamsUtils.memToString("Initial_position", createParams.initialPosition().bytes()));

    final HashBitHandle bitVector = createParams.bitVector();
    if (bitVector != null) {
      sb.append(ParamsUtils.memToString("Bit_vector", bitVector.bytes()));
    }

    long pBytes = 0;
    pBytes += createParams.hash().bytes();
    pBytes += createParams.value().bytes();
    pBytes += createParams.initialPosition().bytes();
    if (bitVector != null) {
      pBytes += bitVector.bytes();
    }
    assert pBytes == bytes(createParams);
  }

  /**
   * Create a human readable description of the memory usage of an <code>IndexImplementation</code>
   * constructed using this.
   * @param createParams parameters that will be used to construct the <code>IndexImplementation</code>.
   * @return the memory description.
   */
  public static String memString(final CreateParams createParams) {
    final StringBuilder sb = new StringBuilder();
    memString(sb, createParams);
    return sb.toString();
  }

  /**
   * Create a human readable description of the memory usage of an <code>IndexImplementation</code>
   * constructed using this.
   * @param sb where to place the result.
   * @param createParams parameters that will be used to construct the <code>IndexImplementation</code>.
   */
  public static void memString(final StringBuilder sb, final CreateParams createParams) {
    sb.append("Memory Usage\tbytes\tlength").append(StringUtils.LS);
    long totalBytes = 0;
    sb.append("\t\t").append(StringUtils.commas(createParams.hash().bytes())).append("\t").append(StringUtils.commas(createParams.hash().length())).append("\tHash").append(StringUtils.LS);
    totalBytes +=  createParams.hash().bytes();

    sb.append("\t\t").append(StringUtils.commas(createParams.value().bytes())).append("\t").append(StringUtils.commas(createParams.value().length())).append("\tValue").append(StringUtils.LS);
    totalBytes +=  createParams.value().bytes();

    sb.append("\t\t").append(StringUtils.commas(createParams.initialPosition().bytes())).append("\t").append(StringUtils.commas(createParams.initialPosition().length())).append("\tInitial Position").append(StringUtils.LS);
    totalBytes += createParams.initialPosition().bytes();

    final HashBitHandle bitVector = createParams.bitVector();
    if (bitVector != null) {
      sb.append("\t\t").append(StringUtils.commas(bitVector.bytes())).append("\t").append(StringUtils.commas(bitVector.length())).append("\tBit vector").append(StringUtils.LS);
      totalBytes += bitVector.bytes();
    }

    sb.append("\t\t").append(StringUtils.commas(totalBytes)).append("\t\tTotal bytes").append(StringUtils.LS);
    assert totalBytes == bytes(createParams);
  }
}

