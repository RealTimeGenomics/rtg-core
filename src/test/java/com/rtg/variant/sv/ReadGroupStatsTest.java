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
    final ReadGroupStats rg1 = new ReadGroupStats("id", 1000);
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
    final String line = rg1.countsString();
    //System.out.println(line);
    final ReadGroupStats rg = new ReadGroupStats(null, line);
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
    final File f = File.createTempFile("rgstats", "txt");
    try {
      FileUtils.stringToFile("id1\t1000\t2\t22\t242\t2\t200\t25000\t2\t200\t25000\t4\t2\t3\t4\n"
                             + "id2\t1000\t2\t22\t242\t2\t200\t20000\t2\t200\t25000\t4\t2\t3\t4\n"
                             + "id2\t1000\t2\t22\t242\t2\t200\t25000\t2\t200\t25000\t4\t2\t3\t4\n"
                             , f);
      Map<String, ReadGroupStats> stats = ReadGroupStats.loadReadGroupStats(null, f);
      assertEquals(2, stats.size());
      assertNotNull(stats.get("id1"));
      assertNotNull(stats.get("id2"));

      final Map<String, String> remap = new HashMap<>();
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
