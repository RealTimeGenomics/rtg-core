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
    this(contigOverlap, Collections.<String, String>emptyMap(), Collections.<String, String>emptyMap());
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
    if (attribute.equals(K_MER_FREQ)) {
      return Integer.toString(kmerFreq(contigId));
    }
    return super.contigAttribute(contigId, attribute);
  }

  @Override
  public void setContigAttribute(long contigId, String attribute, String value) throws IllegalArgumentException, IllegalStateException {
    if (attribute.equals(K_MER_FREQ)) {
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
