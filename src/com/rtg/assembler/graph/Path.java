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
 * The contents of a path.
 */
public interface Path {

  /**
   * @return the number of contigs in the path (&ge; 2).
   */
  int length();

  /**
   * Get the directional contig in the path selected by index.
   * @param index selects the contig in the path.
   * @return a directional contig identifier.
   * The direction (sign) is determined by the direction of the path and the direction of the contig within the path.
   * @throws RuntimeException if  <code>index &lt; 0 or &ge; pathLength(pathId)</code>.
   */
  long contig(int index);

  /**
   * Get the index of the contig within this path.
   * @param contigId directional identifier of contig.
   * @return the position of the contig within the path or -1 if the contig (including its direction) is not in the path.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is &le; 0 or &gt; <code>numberContigs()</code>.
   */
  int index(long contigId) throws IllegalArgumentException;
}
