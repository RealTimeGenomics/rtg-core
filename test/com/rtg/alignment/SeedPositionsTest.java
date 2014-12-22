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
    for (int i = 0; i < seeds.length; i++) {
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
