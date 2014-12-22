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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Contig;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.dna.DNARange;

/**
 */
@TestClass("com.rtg.assembler.graph.implementation.GraphImplementationTest")
final class ContigReverse extends IntegralAbstract implements Contig {
  private final GraphImplementation mGraphImplementation;
  private final long mContig;
  private final long mAcontig;
  private final long mStart;

  /**
   * @param contig the contig identifier.
   * @param graphImplementation the graph that is parent of this (was an inner class before being moved out)
   */
  ContigReverse(GraphImplementation graphImplementation, long contig) {
    mGraphImplementation = graphImplementation;
    mContig = contig;
    mAcontig = -contig;
    mStart = mGraphImplementation.mContigs.get(mAcontig) - 1;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mContig < 0);
    Exam.assertEquals(mContig, -mAcontig);
    return true;
  }

  @Override
  public int length() {
    final long length = mStart - mGraphImplementation.mContigs.get(mAcontig - 1) + 1;
    assert length <= Integer.MAX_VALUE;
    return (int) length;
  }

  @Override
  public byte nt(int index) {
    if (index < 0 || index >= length()) {
      throw new RuntimeException();
    }
    final long nt = mGraphImplementation.mNt.get(mStart - index);
    assert DNARange.RANGE.valid(nt);
    return DNARange.complement((byte) nt);
  }
}
