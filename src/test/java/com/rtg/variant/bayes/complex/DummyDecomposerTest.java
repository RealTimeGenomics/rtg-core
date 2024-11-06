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
