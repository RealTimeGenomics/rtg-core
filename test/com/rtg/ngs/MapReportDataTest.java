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

package com.rtg.ngs;

import java.io.ByteArrayInputStream;
import java.io.IOException;

import com.rtg.ngs.MapReportData.DistributionType;
import com.rtg.sam.SamBamConstants;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.MemoryPrintStream;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class MapReportDataTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testEmptyReporter() {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final MapReportData mr = merger.createMapReportData();
    assertEquals("", mr.toString());

    for (int i = 0; i < 10; i++) {
      final MapReportData mr2 = merger.createMapReportData();
      assertEquals("", mr2.toString());
    }

    assertEquals("", mr.toString());
    final MapReportData mr3 = merger.blendReportData();
    assertEquals("", mr3.toString());
    assertFalse(mr3.isPairedEnd());

    assertEquals(0, mr.getCommandLines().size());
    assertEquals(0, mr3.getCommandLines().size());
  }

  public void testSimpleReporter() {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final MapReportData mr = merger.createMapReportData();
    for (int i = 0; i < 10; i++) {
      final SAMRecord sr = new SAMRecord(null);
      sr.setFlags(0);
      sr.setAttribute("AS", i / 2);
      sr.setAttribute("NM", i / 2);
      sr.setAttribute("NH", i / 2);
      mr.processRead(sr);
    }
    final StringBuilder expected = new StringBuilder();
    expected.append("MAPC\t0\t10").append(StringUtils.LS);
    for (int i = 0; i < 5; i++) {
      expected.append("AS\t").append(i).append("\t2").append(StringUtils.LS);
    }
    expected.append("RLEN\t0\t10").append(StringUtils.LS);
    for (int i = 0; i < 5; i++) {
      expected.append("NH\t").append(i).append("\t2").append(StringUtils.LS);
    }
    for (int i = 0; i < 5; i++) {
      expected.append("NM\t").append(i).append("\t2").append(StringUtils.LS);
    }
    expected.append("ORI\t0\t10").append(StringUtils.LS);
    expected.append("MAPQ\t0\t10").append(StringUtils.LS);

    assertEquals(expected.toString(), mr.toString());

    for (int i = 0; i < 10; i++) {
      final MapReportData mr2 = merger.createMapReportData();
      assertEquals("", mr2.toString());
    }

    assertEquals(expected.toString(), mr.toString());
    final MapReportData mr3 = merger.blendReportData();
    assertEquals(expected.toString(), mr3.toString());
    for (DistributionType type : DistributionType.values()) {
      assertNotNull(mr3.getHistogram(type));
    }
    assertEquals(1, mr3.getHistogram(DistributionType.RLEN).getLength());
    assertFalse(mr3.isPairedEnd());
  }

  public void testSimpleReporter2() {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final MapReportData mr = merger.createMapReportData();
    for (int i = 0; i < 10; i++) {
      final SAMRecord sr = new SAMRecord(null);
      sr.setFlags(SamBamConstants.SAM_READ_IS_PAIRED);
      sr.setAttribute("AS", i / 2);
      sr.setAttribute("NM", i / 2);
      sr.setAttribute("NH", i / 2);
      sr.setMappingQuality(20 / (i % 2 + 1));
      mr.processRead(sr);
    }
    final StringBuilder expected = new StringBuilder();
    expected.append("MAPC\t0\t10").append(StringUtils.LS);
    expected.append("MAPC\t2\t10").append(StringUtils.LS);
    for (int i = 0; i < 5; i++) {
      expected.append("AS\t").append(i).append("\t2").append(StringUtils.LS);
    }
    expected.append("RLEN\t0\t10").append(StringUtils.LS);
    for (int i = 0; i < 5; i++) {
      expected.append("NH\t").append(i).append("\t2").append(StringUtils.LS);
    }
    for (int i = 0; i < 5; i++) {
      expected.append("NM\t").append(i).append("\t2").append(StringUtils.LS);
    }
    expected.append("ORI\t0\t10").append(StringUtils.LS);
    expected.append("MAPQ\t10\t5").append(StringUtils.LS);
    expected.append("MAPQ\t20\t5").append(StringUtils.LS);

    assertEquals(expected.toString(), mr.toString());

    for (int i = 0; i < 10; i++) {
      final MapReportData mr2 = merger.createMapReportData();
      assertEquals("", mr2.toString());
    }

    assertEquals(expected.toString(), mr.toString());
    final MapReportData mr3 = merger.blendReportData();
    assertEquals(expected.toString(), mr3.toString());
    assertFalse(mr3.isPairedEnd());
  }

  public void testUnmappedReads() {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final MapReportData mr = merger.createMapReportData();
    for (int i = 0; i < 10; i++) {
      final SAMRecord sr = new SAMRecord(null);
      sr.setFlags(SamBamConstants.SAM_READ_IS_UNMAPPED);
      sr.setAttribute("AS", i / 2);
      sr.setAttribute("NM", i / 2);
      sr.setAttribute("NH", i / 2);
      mr.processRead(sr);
    }
    final StringBuilder expected = new StringBuilder();
    expected.append("MAPC\t3\t10").append(StringUtils.LS);
    expected.append("RLENU\t0\t10").append(StringUtils.LS);

    assertEquals(expected.toString(), mr.toString());

    for (int i = 0; i < 10; i++) {
      final MapReportData mr2 = merger.createMapReportData();
      assertEquals("", mr2.toString());
    }

    assertEquals(expected.toString(), mr.toString());
    final MapReportData mr3 = merger.blendReportData();
    assertEquals(expected.toString(), mr3.toString());
    assertFalse(mr3.isPairedEnd());
  }

  public void testPairedReads() {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final MapReportData mr = merger.createMapReportData();
    for (int i = 0; i < 4; i++) {
      final SAMRecord sr = new SAMRecord(null);
      int flag = (i % 2 == 0) ? SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR : SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR;
      if (i % 2 == 0) {
        flag |= SamBamConstants.SAM_READ_IS_REVERSE;
      } else {
        flag |= SamBamConstants.SAM_MATE_IS_REVERSE;
      }
      if (i > 2) {
        flag |= SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR;
        sr.setInferredInsertSize((i % 2 == 0) ? -123 : 123);
      }
      flag |= SamBamConstants.SAM_READ_IS_PAIRED;
      sr.setFlags(flag);
      sr.setAttribute("AS", i / 2);
      sr.setAttribute("NM", i / 2);
      sr.setAttribute("NH", i / 2);
      mr.processRead(sr);
    }
    final StringBuilder expected = new StringBuilder();
    expected.append("MAPC\t0\t4").append(StringUtils.LS);
    expected.append("MAPC\t1\t1").append(StringUtils.LS);
    expected.append("MAPC\t2\t3").append(StringUtils.LS);
    for (final String s : new String[]{"AS\t0\t1", "AS\t1\t1", "AS2\t0\t1", "AS2M\t1\t1", "RLEN\t0\t2", "RLEN2\t0\t1", "RLEN2M\t0\t1", "NH\t0\t1", "NH\t1\t1", "NH2\t0\t1", "NH2M\t1\t1", "NM\t0\t1", "NM\t1\t1", "NM2\t0\t1", "NM2M\t1\t1", "ORI\t1\t2", "ORI2\t0\t1", "ORI2M\t0\t1", "MORI\t1\t1", "FLEN\t123\t1"}) {
      expected.append(s).append(StringUtils.LS);
    }
    expected.append("MAPQ\t0\t4").append(StringUtils.LS);

    assertEquals(expected.toString(), mr.toString());

    for (int i = 0; i < 10; i++) {
      final MapReportData mr2 = merger.createMapReportData();
      assertEquals("", mr2.toString());
    }

    assertEquals(expected.toString(), mr.toString());
    final MapReportData mr3 = merger.blendReportData();
    assertEquals(expected.toString(), mr3.toString());
    assertTrue(mr3.isPairedEnd());
  }


  public void testMultipleReporters() {
    final MapReportData.Merger merger = new MapReportData.Merger();
    for (int i = 0; i < 10; i++) {
      final MapReportData mr = merger.createMapReportData();
      final SAMRecord sr = new SAMRecord(null);
      sr.setFlags(SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR | SamBamConstants.SAM_READ_IS_PAIRED);
      sr.setAttribute("AS", i / 2);
      sr.setInferredInsertSize((i % 3 == 0 ? -1 : 1) * (i + 42));
      sr.setReadString("ACGT");
      mr.processRead(sr);
      final String str = mr.toString();
      //System.err.println(str);
      assertTrue(str.contains("ASM\t" + (i / 2) + "\t1" + StringUtils.LS));
      assertTrue(str.contains("RLENM\t4\t1" + StringUtils.LS));
    }
    final StringBuilder expected = new StringBuilder();
    expected.append("MAPC\t0\t10").append(StringUtils.LS);
    expected.append("MAPC\t1\t10").append(StringUtils.LS);
    for (int i = 0; i < 5; i++) {
      expected.append("ASM\t").append(i).append("\t2").append(StringUtils.LS);
    }
    expected.append("RLENM\t4\t10").append(StringUtils.LS);
    expected.append("ORIM\t0\t10").append(StringUtils.LS);
    expected.append("MORI\t0\t6").append(StringUtils.LS);
    for (int i = 0; i < 10; i++) {
      if ((i % 3) != 0) {
        expected.append("FLEN\t").append(i + 42).append("\t1").append(StringUtils.LS);
      }
    }
    expected.append("MAPQ\t0\t10").append(StringUtils.LS);

    final MapReportData mr3 = merger.blendReportData();
    assertEquals(expected.toString(), mr3.toString());
    assertTrue(mr3.isPairedEnd());
  }

  public void testReadWrite() throws IOException {
    final MapReportData.Merger merger = new MapReportData.Merger();
    createMapData(merger);

    final MapReportData mr3 = merger.blendReportData();
    final MemoryPrintStream mps = new MemoryPrintStream();
    assertEquals(0, mr3.getCommandLines().size());
    mr3.write(mps.printStream());

    final String out = mps.toString();
    assertNotNull(out);
    assertTrue(out.contains("#Version\t"));
    assertTrue(out.contains(MapReportData.VERSION));
    TestUtils.containsAll(out, "#Alignment Score Distribution", "#Fragment Length Distribution", "#Mapping Quality (MAPQ) Distribution");
    mps.printStream().println("#CL\tfoooooo");

    final MapReportData mr4 = merger.createMapReportData(new ByteArrayInputStream(mps.outputStream().toByteArray()));
    assertEquals(mr3.toString(), mr4.toString());
    assertFalse(mr3.isPairedEnd());
    assertFalse(mr4.isPairedEnd());
    assertEquals(2, mr4.getCommandLines().size());
    assertTrue(mr4.getCommandLines().contains("foooooo"));
  }

  public static void createMapData(MapReportData.Merger merger) {
    for (int i = 0; i < 10; i++) {
      final MapReportData mr = merger.createMapReportData();
      final SAMRecord sr = new SAMRecord(null);
      sr.setFlags(SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR);
      sr.setAttribute("AS", i / 2);
      sr.setInferredInsertSize(i / 2 + 42);
      mr.processRead(sr);
    }
  }

  public void testReadWriteEmpty() throws IOException {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final MapReportData mr3 = merger.blendReportData();
    final MemoryPrintStream mps = new MemoryPrintStream();
    mr3.write(mps.printStream());

    final String out = mps.toString();
    assertNotNull(out);
    final String[] lines = out.split(StringUtils.LS);
    assertEquals(28, lines.length);
    for (String line : lines) {
      assertTrue(line.startsWith("#"));
    }
    //assertEquals("#VERSION: 0" + StringUtils.LS + "#AS distribution" + StringUtils.LS + "#FLEN distribution" + StringUtils.LS, out);

    final MapReportData mr4 = merger.createMapReportData(new ByteArrayInputStream(mps.outputStream().toByteArray()));
    assertEquals(mr3.toString(), mr4.toString());
    assertFalse(mr3.isPairedEnd());
    assertFalse(mr4.isPairedEnd());
  }

  public void testBadVersion() throws IOException {
    final MapReportData.Merger merger = new MapReportData.Merger();
    final ByteArrayInputStream bais = new ByteArrayInputStream("#Version: 3\n".getBytes());
    try {
      merger.createMapReportData(bais);
      fail("Expected Exception");
    } catch (final NoTalkbackSlimException e) {
      assertTrue(e.getMessage().contains("Unsupported map statistics version")); // expected
    }
  }

  public void testDistributionType() {
    TestUtils.testEnum(DistributionType.class, "[MAPC, AS, AS2, ASM, AS2M, RLEN, RLEN2, RLENM, RLEN2M, RLENU, NH, NH2, NHM, NH2M, NM, NM2, NMM, NM2M, ORI, ORI2, ORIM, ORI2M, MORI, FLEN, MAPQ]");
    assertEquals(DistributionType.AS, DistributionType.AS.getType(false, false));
    assertEquals(DistributionType.AS2, DistributionType.AS.getType(true, false));
    assertEquals(DistributionType.ASM, DistributionType.AS.getType(false, true));
    assertEquals(DistributionType.AS2M, DistributionType.AS.getType(true, true));

    for (final DistributionType dt : DistributionType.values()) {
      final String ln = dt.longName();
      assertTrue(ln.contains("Score") || ln.contains("Length") || ln.contains("Orientation") || ln.contains("Hits") || ln.contains("Mismatches") || ln.contains("Quality") || ln.contains("Counts"));
    }
  }
}
