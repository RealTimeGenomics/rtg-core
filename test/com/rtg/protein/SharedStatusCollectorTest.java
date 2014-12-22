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
package com.rtg.protein;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;

import com.rtg.ngs.MapStatisticsArm;
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

    ByteArrayOutputStream bos = new ByteArrayOutputStream();
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
    assertEquals(22, stats.value(MapStatisticsField.TOTAL_READS, MapStatisticsArm.LEFT));
    assertEquals(1, stats.value(MapStatisticsField.UNMATED_UNIQUE_READS, MapStatisticsArm.LEFT));
    assertEquals(20, stats.value(MapStatisticsField.UNMAPPED_UNMATED_POOR, MapStatisticsArm.LEFT));
    assertEquals(1, stats.value(MapStatisticsField.UNMAPPED_NO_HITS, MapStatisticsArm.LEFT));
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
            for (int j = 0; j < 100; j++) {
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
      for (int i = 0; i < 100; i++) {
        assertEquals(0x3f, ssc.getStatus(i));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

}
