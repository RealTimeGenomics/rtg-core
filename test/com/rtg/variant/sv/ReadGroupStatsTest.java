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

package com.rtg.variant.sv;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import com.rtg.util.io.FileUtils;

import junit.framework.TestCase;

/**
 */
public class ReadGroupStatsTest extends TestCase {


  public void testParseSingle() {
    try {
      new ReadGroupStats(null, "");
      fail();
    } catch (IllegalArgumentException e) {
    }
    try {
      new ReadGroupStats(null, "id\t1\t2\t3\t4\t5foo\t7\t8\t9");
      fail();
    } catch (IllegalArgumentException e) {
    }
    ReadGroupStats rg1 = new ReadGroupStats("id", 1000);
    rg1.addLength(11);
    rg1.addLength(11);
    rg1.addProper();
    rg1.addGapSize(50);
    rg1.addProper();
    rg1.addGapSize(150);
    rg1.addScore(2);
    rg1.addScore(4);
    rg1.addDiscordant();
    rg1.addDiscordant();
    rg1.addDiscordant();
    rg1.addUnmated();
    rg1.addUnmated();
    rg1.addUnmated();
    rg1.addUnmated();
    String line = rg1.countsString();
    //System.out.println(line);
    ReadGroupStats rg = new ReadGroupStats(null, line);
    //System.out.println(rg.toString());

    assertEquals("id", rg.id());
    assertEquals(11.0, rg.meanLength());
    assertEquals(100.0, rg.gapMean());
    assertEquals(70.71, rg.gapStdDev(), 0.01);
    assertEquals(4, rg.maxAlignment());
    assertEquals(0.0015, rg.properRate(), 0.0001);
    assertEquals(0.000005, rg.properRandomRate(), 0.000001);
    assertEquals(0.002, rg.discordantRate(), 0.001);
    assertEquals(0.002, rg.unmatedRate(), 0.001);
  }
  public void testParseFile() throws IOException {
    File f = File.createTempFile("rgstats", "txt");
    try {
      FileUtils.stringToFile("id1\t1000\t2\t22\t242\t2\t200\t25000\t2\t200\t25000\t4\t2\t3\t4\n"
                             + "id2\t1000\t2\t22\t242\t2\t200\t20000\t2\t200\t25000\t4\t2\t3\t4\n"
                             + "id2\t1000\t2\t22\t242\t2\t200\t25000\t2\t200\t25000\t4\t2\t3\t4\n"
                             , f);
      Map<String, ReadGroupStats> stats = ReadGroupStats.loadReadGroupStats(null, f);
      assertEquals(2, stats.size());
      assertNotNull(stats.get("id1"));
      assertNotNull(stats.get("id2"));

      Map<String, String> remap = new HashMap<>();
      remap.put("id1", "merged");
      remap.put("id2", "merged");
      stats = ReadGroupStats.loadReadGroupStats(remap, f);
      assertEquals(1, stats.size());
      assertNotNull(stats.get("merged"));
      assertNull(stats.get("id2"));
    } finally {
      assertTrue(f.delete());
    }

  }
}
