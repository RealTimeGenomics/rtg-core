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
package com.rtg.protein;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;

import com.rtg.reader.Arm;
import com.rtg.ngs.MapStatisticsField;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class SharedStatusCollectorTest extends TestCase {

  /**
   */
  public SharedStatusCollectorTest(String name) {
    super(name);
  }

  public void testUnmatedWriting() throws Exception {
    final SharedStatusCollector res = new SharedStatusCollector(1, null);
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    res.writeUnmapped(bos, null, 0);
    assertEquals("0" + LS, bos.toString());

    bos.reset();
    res.setStatus(0, SharedStatusCollector.EXCEEDS_ALIGNMENT_THRESHOLD);
    res.writeUnmapped(bos, null, 0);
    assertEquals("0\td" + LS, bos.toString());

    bos.reset();
    res.setStatus(0, SharedStatusCollector.EXCEEDS_PERCENT_ID_THRESHOLD);
    res.writeUnmapped(bos, null, 0);
    assertEquals("0\tf" + LS, bos.toString());

    bos.reset();
    res.setStatus(0, SharedStatusCollector.EXCEEDS_E_SCORE_THRESHOLD);
    res.writeUnmapped(bos, null, 0);
    assertEquals("0\tg" + LS, bos.toString());

    bos.reset();
    res.setStatus(0, SharedStatusCollector.BELOW_BIT_SCORE_THRESHOLD);
    res.writeUnmapped(bos, null, 0);
    assertEquals("0\th" + LS, bos.toString());

    bos.reset();
    res.setStatus(0, SharedStatusCollector.EXCEEDS_N_THRESHOLD);
    res.writeUnmapped(bos, null, 0);
    assertEquals("0\te" + LS, bos.toString());

    bos.reset();
    res.setStatus(0, SharedStatusCollector.RESULT_WRITTEN);
    res.writeUnmapped(bos, null, 0);
    assertEquals("", bos.toString());
  }

  public void testCalculateStatistics() {
    final MapXStatistics stats = new MapXStatistics(null);
    final SharedStatusCollector collector = new SharedStatusCollector(22, stats);
    collector.setStatus(0, (byte) 0x20);
    collector.setStatus(1, (byte) 0x01);
    collector.setStatus(2, (byte) 0x01);
    collector.setStatus(3, (byte) 0x02);
    collector.setStatus(4, (byte) 0x02);
    collector.setStatus(5, (byte) 0x02);
    collector.setStatus(6, (byte) 0x04);
    collector.setStatus(7, (byte) 0x04);
    collector.setStatus(8, (byte) 0x04);
    collector.setStatus(9, (byte) 0x04);
    collector.setStatus(10, (byte) 0x08);
    collector.setStatus(11, (byte) 0x08);
    collector.setStatus(12, (byte) 0x08);
    collector.setStatus(13, (byte) 0x08);
    collector.setStatus(14, (byte) 0x08);
    collector.setStatus(15, (byte) 0x10);
    collector.setStatus(16, (byte) 0x10);
    collector.setStatus(17, (byte) 0x10);
    collector.setStatus(18, (byte) 0x10);
    collector.setStatus(19, (byte) 0x10);
    collector.setStatus(20, (byte) 0x10);
    collector.calculateStatistics();
    assertEquals(22, stats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));
    assertEquals(1, stats.value(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(20, stats.value(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT));
    assertEquals(1, stats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
    final byte[] testGetReadProtein = collector.getReadProtein(131);
    assertNull(testGetReadProtein);
  }

  public void testMt() throws Exception {
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    try {
      final MapXStatistics stats = new MapXStatistics(null);
      final SharedStatusCollector ssc = new SharedStatusCollector(100, stats);
      final SimpleThreadPool stp = new SimpleThreadPool(16, "TestSharedStatusMt", true);

      for (byte i = 1; i < 0x21; i = (byte) (i << 1)) {
        final byte status = i;
        stp.execute(new IORunnable() {
          @Override
          public void run() {
            for (int j = 0; j < 100; ++j) {
              ssc.setStatus(j, status);
              try {
              Thread.sleep(1);
              } catch (Exception e) {
              }
            }
          }
        });
      }
      stp.terminate();
      for (int i = 0; i < 100; ++i) {
        assertEquals(0x3f, ssc.getStatus(i));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

}
