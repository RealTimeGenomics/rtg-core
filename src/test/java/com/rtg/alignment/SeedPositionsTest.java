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
package com.rtg.alignment;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import com.rtg.mode.DnaUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class SeedPositionsTest extends TestCase {

  public void testSeedPositions1() {
    final SeedPositions s = new SeedPositions();
    s.mX1 = 1;
    s.mX2 = 2;
    s.mY1 = 3;
    s.mY2 = 5;
    assertEquals("[1 2) [3 5) ", s.toString());
    assertEquals(1, s.xWidth());
    assertEquals(2, s.yWidth());
  }

  private SeedPositions[] getSeeds() {
    final SeedPositions[] seeds = new SeedPositions[3];
    for (int i = 0; i < seeds.length; ++i) {
      seeds[i] = new SeedPositions();
      seeds[i].mType = 0;
      final int start = i * i;
      final int end = (i + 1) * (i + 1);
      seeds[i].mX1 = start;
      seeds[i].mX2 = end;
      seeds[i].mY1 = start + i;
      seeds[i].mY2 = end + i;
    }
    return seeds;
  }

  static final String DUMP_SEEDS1 = ""
    + "[0 1) [0 1)  A A" + StringUtils.LS
    + "[1 4) [2 5)  CCT CCT" + StringUtils.LS
    + "[4 9) [6 11)  AACCG AACCG" + StringUtils.LS;
  static final byte[] READ1 = DnaUtils.encodeString("ACCTAACCG");
  static final byte[] TMPL1 = DnaUtils.encodeString("AcCCTgAACCGG");

  public void testDumpSeeds() {
    final SeedPositions[] seeds = getSeeds();
    assertEquals(""
        + "[0 1) [0 1) " + StringUtils.LS
        + "[1 4) [2 5) " + StringUtils.LS
        + "[4 9) [6 11) " + StringUtils.LS, SeedPositions.dumpSeedsAsString(seeds));
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream orig = System.err;
    System.setErr(new PrintStream(bos));
    try {
      SeedPositions.dumpSeeds(seeds, seeds.length, READ1, TMPL1);
      System.err.flush();
      assertEquals(DUMP_SEEDS1, bos.toString());
    } finally {
      System.setErr(orig);
    }
  }

  public void testSeedIntegrity() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream orig = System.err;
    System.setErr(new PrintStream(bos)); // ignore stderr output
    try {
      SeedPositions[] seeds = getSeeds();
      SeedPositions.seedIntegrity(seeds, seeds.length, READ1, READ1.length, TMPL1, 0);

      seeds = getSeeds();
      seeds[0].mX2 = seeds[0].mX1 - 1;
      checkException(seeds, "x -ve range");

      seeds = getSeeds();
      seeds[1].mY2 = seeds[1].mY1 - 1;
      checkException(seeds, "x and y ranges are different lengths");

      seeds = getSeeds();
      seeds[0].mX2 = seeds[0].mX1;
      seeds[0].mY2 = seeds[0].mY1;
      checkException(seeds, "empty seed: [0 0) [0 0) ");

      seeds = getSeeds();
      seeds[1].mX1--;
      seeds[1].mY1--;
      checkException(seeds, "non monotonic x: [0 1) [0 1)  [0 4) [1 5)  9 0");

      // several off-end errors
      seeds = getSeeds();
      final String offEndMsg = "BUG: off end of array";
      seeds[0].mX1--;
      seeds[0].mX2--;
      checkException(seeds, offEndMsg);
      seeds = getSeeds();
      seeds[0].mY1--;
      seeds[0].mY2--;
      checkException(seeds, offEndMsg);
      seeds = getSeeds();
      seeds[seeds.length - 1].mX1++;
      seeds[seeds.length - 1].mX2++;
      checkException(seeds, offEndMsg);
      seeds = getSeeds();
      seeds[seeds.length - 1].mY1 += 2;
      seeds[seeds.length - 1].mY2 += 2;
      checkException(seeds, offEndMsg);

      // unequal contents messages in developer's log
      final ByteArrayOutputStream logbos = new ByteArrayOutputStream();
      final PrintStream pos = new PrintStream(logbos);
      Diagnostic.setLogStream(pos);
      seeds = getSeeds();
      seeds[0].mY1++;
      seeds[0].mY2++;
      SeedPositions.seedIntegrity(seeds, seeds.length, READ1, READ1.length, TMPL1, 0);
      seeds = getSeeds();
      seeds[1].mY1--;
      seeds[1].mY2--;
      SeedPositions.seedIntegrity(seeds, seeds.length, READ1, READ1.length, TMPL1, 0);
      Diagnostic.setLogStream();
      final String log = logbos.toString();
      //System.out.println(log);
      assertTrue(log.contains("bad seed at 0 1   values = 1 2 [0 1) [1 2)  "));
      assertTrue(log.contains(StringUtils.LS
          + "ACCTAACCG" + StringUtils.LS
          + "ACCCTGAAC" + StringUtils.LS
          + "A" + StringUtils.LS
          + "C" + StringUtils.LS));
      assertTrue(log.contains("bad seed at 3 3   values = 4 2 [1 4) [1 4)  "));
      assertTrue(log.contains(StringUtils.LS
          + "ACCTAACCG" + StringUtils.LS
          + "ACCCTGAAC" + StringUtils.LS
          + "CCT" + StringUtils.LS
          + "CCC" + StringUtils.LS));
    } finally {
      System.setErr(orig);
    }
  }

  public void testOverlapIntegrity() {
    final SeedPositions[] seeds = getSeeds();
    assertTrue(SeedPositions.overlapIntegrity(seeds, seeds.length, READ1, READ1.length, TMPL1, 0));
    seeds[1].mX1--;
    seeds[1].mY1--;
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());
      try {
        assertFalse(SeedPositions.overlapIntegrity(seeds, seeds.length, READ1, READ1.length, TMPL1, 0));
        assertTrue(mps.toString(), mps.toString().contains("overlapping seeds"));
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

  private void checkException(final SeedPositions[] seeds, final String msg) {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());
      try {
        assertFalse(SeedPositions.seedIntegrity(seeds, seeds.length, READ1, READ1.length, TMPL1, 0));
        assertTrue(mps.toString(), mps.toString().contains(msg));
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }
}
