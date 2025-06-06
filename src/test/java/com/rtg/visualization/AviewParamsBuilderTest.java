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

package com.rtg.visualization;

import java.io.File;

import junit.framework.TestCase;

/**
 */
public class AviewParamsBuilderTest extends TestCase {

  /**
   * Test method for {@link com.rtg.visualization.AviewParamsBuilder#create()}.
   */
  public final void testCreateDefault() {
    final AviewParams p = new AviewParamsBuilder().create();
    assertEquals(true, p.displayDots());
    assertEquals(true, p.useTerminalColor());
    assertEquals(false, p.printCigars());
    assertEquals(false, p.printReadName());

    assertEquals(-1, p.start());
    assertEquals(-1, p.end());
    assertEquals(-1, p.regionPadding());
    assertEquals(5, p.headerLineRepeat());
    assertEquals(Integer.MAX_VALUE, p.maxMatedAlignmentScore());
    assertEquals(Integer.MAX_VALUE, p.maxUnmatedAlignmentScore());
    assertEquals(Integer.MAX_VALUE, p.maxIhScore());

    assertEquals("", p.sequenceName());

    assertEquals(0, p.alignmentsFiles().length);
    assertEquals(0, p.unmappedFiles().length);
    assertEquals(0, p.readFiles().length);
    assertEquals(0, p.trackFiles().length);
    assertNull(p.baselineFile());
    assertNull(p.referenceFile());
    assertFalse(p.sortReads());
  }

  /**
   * Test method for {@link com.rtg.visualization.AviewParamsBuilder#create()}.
   */
  public final void testCreateValues() {
    final AviewParams p = new AviewParamsBuilder()
    .displayDots(false)
    .useTerminalColor(false)
    .printCigars(true)
    .printReadName(true)
    .region("foo:10-20")
    .headerLineRepeat(6)
    .maxIhScore(4)
    .maxMatedAlignmentScore(1)
    .maxUnmatedAlignmentScore(2)
    .alignments(new File[] {new File("a"), new File("b")})
    .unmapped(new File[] {new File("una"), new File("unb")})
    .reads(new File[] {new File("ra"), new File("rb")})
    .trackFiles(new File[]{new File("foundSnps")})
    .baselineFile(new File("generatedSnps"))
    .reference(new File("template"))
    .regionPadding(40)
    .sortReads(true)
    .create();
    assertEquals(false, p.displayDots());
    assertEquals(false, p.useTerminalColor());
    assertEquals(true, p.printCigars());
    assertEquals(true, p.printReadName());

    assertEquals(10, p.start());
    assertEquals(21, p.end());
    assertEquals(6, p.headerLineRepeat());
    assertEquals(1, p.maxMatedAlignmentScore());
    assertEquals(2, p.maxUnmatedAlignmentScore());
    assertEquals(4, p.maxIhScore());
    assertEquals(40, p.regionPadding());
    assertEquals("foo", p.sequenceName());

    assertEquals("a", p.alignmentsFiles()[0].getName());
    assertEquals("b", p.alignmentsFiles()[1].getName());

    assertEquals("una", p.unmappedFiles()[0].getName());
    assertEquals("unb", p.unmappedFiles()[1].getName());

    assertEquals("ra", p.readFiles()[0].getName());
    assertEquals("rb", p.readFiles()[1].getName());

    assertEquals("foundSnps", p.trackFiles()[0].getName());
    assertEquals("generatedSnps", p.baselineFile().getName());
    assertEquals("template", p.referenceFile().getName());

    assertTrue(p.sortReads());
  }

}

