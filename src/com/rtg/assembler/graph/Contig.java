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

package com.rtg.assembler.graph;

/**
 * The contents of a contig.
 */
public interface Contig {

  /**
   * @return the number of nucleotides in the contig.
   */
  int length();

  /**
   * Get the specified nucleotide.
   * @param index selects the nucleotide.
   * @return the specified nucleotide.
   * @throws RuntimeException if <code>index &lt; 0 or &ge; length()</code>.
   */
  byte nt(int index);
}
