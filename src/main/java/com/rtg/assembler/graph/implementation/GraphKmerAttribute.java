/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.assembler.graph.implementation;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.integrity.Exam;

/**
 */
public class GraphKmerAttribute extends GraphImplementation {

  /** Attribute key for the <code>k-mer</code> frequency. */
  public static final String K_MER_FREQ = "kMerFreq";
  /** Attribute for representing the number of reads mapping to a contig or path */
  public static final String READ_COUNT = "readCount";
  /** Attribute description for read count */
  public static final String READ_COUNT_DESCRIPTION = "number of reads that map here";

  private static Map<String, String> makeContigAttributes(Map<String, String> contigAttributes) {
    final Map<String, String> contigA = new HashMap<>(contigAttributes);
    contigA.put(K_MER_FREQ, "Sum of frequencies of the k-mers in the contig");
    return contigA;
  }

  private final ExtensibleIndex mKmerFreq = new IntChunks(1);

  /**
   * Construct with no additional attributes.
   * @param contigOverlap the number of bases adjacent contigs overlap
   */
  public GraphKmerAttribute(int contigOverlap) {
    this(contigOverlap, Collections.emptyMap(), Collections.emptyMap());
  }

  /**
   * @param contigOverlap the number of bases adjacent contigs overlap
   * @param contigAttributes keys and comments for contig attributes.
   * @param pathAttributes keys for path attributes.
   */
  public GraphKmerAttribute(int contigOverlap, Map<String, String> contigAttributes, Map<String, String> pathAttributes) {
    super(contigOverlap, makeContigAttributes(contigAttributes), pathAttributes);
  }

  /**
   * @param contigId signed contig identifier.
   * @return the value of the <code>kMerFreq</code>.
   */
  public int kmerFreq(final long contigId) {
    final long acontig = absContig(contigId);
    return (int) mKmerFreq.get(acontig);
  }

  /**
   * @param contigId signed contig identifier.
   * @param value to set the attribute to.
   */
  public void setKmerFreq(final long contigId, final int value) {
    assert 0 <= value;
    final long acontig = absContig(contigId);
    mKmerFreq.set(acontig, value);
  }

  @Override
  protected void contigExpand() {
    super.contigExpand();
    mKmerFreq.extendBy(1);
  }

  @Override
  public String contigAttribute(long contigId, String attribute) {
    if (K_MER_FREQ.equals(attribute)) {
      return Integer.toString(kmerFreq(contigId));
    }
    return super.contigAttribute(contigId, attribute);
  }

  @Override
  public void setContigAttribute(long contigId, String attribute, String value) {
    if (K_MER_FREQ.equals(attribute)) {
      final int ival = Integer.parseInt(value);
      setKmerFreq(contigId, ival);
    }
    super.setContigAttribute(contigId, attribute, value);
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertEquals(numberContigs(), mKmerFreq.length() - 1);
    return true;
  }
}
