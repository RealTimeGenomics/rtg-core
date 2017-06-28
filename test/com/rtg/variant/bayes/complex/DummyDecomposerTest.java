/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.bayes.complex;

import java.util.Arrays;

import com.rtg.variant.VariantLocus;

import junit.framework.TestCase;

/**
 */
public class DummyDecomposerTest extends TestCase {

  public void testAnchorBase() {
    final VariantLocus locus = new VariantLocus("blah", 0, 3, "CGT", 'A');
    assertEquals('A', AbstractDecomposer.getAnchorBase(locus, 0));
    assertEquals('C', AbstractDecomposer.getAnchorBase(locus, 1));
    assertEquals('G', AbstractDecomposer.getAnchorBase(locus, 2));
    assertEquals('T', AbstractDecomposer.getAnchorBase(locus, 3));
  }

  public void testExtractAlts() {
    assertEquals("[CCCC]", AbstractDecomposer.extractAlts(AbstractDecomposerTest.createVariant("ACGTCTGTCT", "CCCC", null)).toString());
    assertEquals("[CCCC, AC]", AbstractDecomposer.extractAlts(AbstractDecomposerTest.createVariant("ACGTCTGTCT", "CCCC", "AC")).toString());
  }

  public void testExtractAlleles() {
    assertEquals("[ACGTCTGTCT, CCCC]", Arrays.toString(AbstractDecomposer.extractAlleles(AbstractDecomposerTest.createVariant("ACGTCTGTCT", "CCCC", null))));
    assertEquals("[ACGTCTGTCT, CCCC, AC]", Arrays.toString(AbstractDecomposer.extractAlleles(AbstractDecomposerTest.createVariant("ACGTCTGTCT", "CCCC", "AC"))));
  }
}
