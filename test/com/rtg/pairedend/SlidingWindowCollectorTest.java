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
package com.rtg.pairedend;


import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Properties;
import java.util.Random;

import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.tempstage.PairedTempFileWriter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class SlidingWindowCollectorTest extends TestCase {

  private File mDir;

  @Override
  public void setUp() throws Exception {
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    FileHelper.deleteAll(mDir);
    mDir = null;
  }

  private static class TestPairedAlignmentWriter implements PairedTempFileWriter {
    private final int mExpectedTemplateCount;
    private int mTemplateCount = 0;

    private final int mExpectedPairResultCount;
    private int mPairResultCount = 0;
    private int mLeftPosition; // initialised by nextTemplateId

    private int mExpectedIndex = 0;
    private long[] mExpectedReadIds = null;
    private boolean[] mExpectedFirsts = null;
    private boolean[] mExpectedRC1s = null;
    private int[] mExpectedStart1s = null;

    private int mPairedLeft = 0;

    // Implementation of com.rtg.pairedend.PairedAlignmentWriter

    @Override
    public void initialiseUnmated(OutputStream unmappedOut, MapQScoringReadBlocker a, MapQScoringReadBlocker b) {
    }
    @Override
    public void initialiseMated(OutputStream unmappedOut) {
    }
    @Override
    public void closeMated() {
    }

    @Override
    public void closeUnmated() {
    }

    TestPairedAlignmentWriter(int expectedTemplateCount) {
      this(expectedTemplateCount, 0);
    }

    @Override
    public boolean unmatedResult(int readId, boolean first, boolean rc, int start) {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    TestPairedAlignmentWriter(int expectedTemplateCount, int expectedPairResultCount) {
      mExpectedTemplateCount = expectedTemplateCount;
      mExpectedPairResultCount = expectedPairResultCount;
    }

    public void setExpectedReadIds(final long[] readIds) {
      mExpectedReadIds = readIds;
    }

    public void setExpectedFirsts(final boolean[] firsts) {
      mExpectedFirsts = firsts;
    }

    public void setExpectedRC1s(final boolean[] rc1s) {
      mExpectedRC1s = rc1s;
    }

    public void setExpectedStart1s(final int[] start1s) {
      mExpectedStart1s = start1s;
    }

    public void pairResult(final int readId, final boolean first, final boolean rc1, final int start1) {
//      System.err.println("PR: " + readId + " : " + first + " : " + rc1 + " : " + start1 + " : " + rc2 + " : " + start2);
      mPairResultCount++;
      // check successive calls produce output in non-descending order
      assertTrue(start1 + " : " + mLeftPosition, start1 >= mLeftPosition);
      mLeftPosition = start1;

      if (mExpectedReadIds != null) {
        assertEquals(mExpectedReadIds[mExpectedIndex], readId);
      }
      if (mExpectedFirsts != null) {
        assertEquals(mExpectedFirsts[mExpectedIndex], first);
      }
      if (mExpectedRC1s != null) {
        assertEquals(mExpectedRC1s[mExpectedIndex], rc1);
      }
      if (mExpectedStart1s != null) {
        assertEquals(mExpectedStart1s[mExpectedIndex], start1);
      }
      mExpectedIndex++;
    }

    @Override
    public boolean pairResultLeft(final MatedHitInfo matedHitInfo) {
      mPairedLeft++;
      pairResult(matedHitInfo.mReadId, !matedHitInfo.mFirstRight, matedHitInfo.mReverseComplementLeft, matedHitInfo.mTemplateStartLeft
      );
      return true;
    }

    @Override
    public void pairResultRight(final MatedHitInfo matedHitInfo) throws IOException {
      pairResult(matedHitInfo.mReadId, matedHitInfo.mFirstRight, matedHitInfo.mReverseComplementRight, matedHitInfo.mTemplateStartRight
      );
    }

    @Override
    public void nextTemplateId(final long templateId) {
      mTemplateCount++;
      mLeftPosition = Integer.MIN_VALUE;
    }

    @Override
    public void close() {
      assertEquals("nextTemplateCount != " + mExpectedTemplateCount, mExpectedTemplateCount, mTemplateCount);
      assertTrue(mPairResultCount % 2 == 0);
      assertEquals("pairResultCount != " + mExpectedPairResultCount, mExpectedPairResultCount, mPairResultCount);
    }
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
  + "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" + StringUtils.LS;
  static final String READ_LEFT_LENGTH2 = ">r1" + StringUtils.LS + "aa" + StringUtils.LS
  + ">r2" + StringUtils.LS + "aa" + StringUtils.LS;
  static final String READ_RIGHT_LENGTH2 = ">r1" + StringUtils.LS + "aa" + StringUtils.LS
  + ">r2" + StringUtils.LS + "aa" + StringUtils.LS;
  static final String AS50 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  //                          12345678901234567890123456789012345678901234567890
  static final String READ_LEFT_LENGTH50 = ">r1" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r2" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r3" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r4" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r5" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r6" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r7" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r8" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r9" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r10" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r11" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r12" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r13" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r14" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r15" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r16" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r17" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r18" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r19" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r20" + StringUtils.LS + AS50 + StringUtils.LS
  + ">r21" + StringUtils.LS + AS50 + StringUtils.LS;
  static final String READ_RIGHT_LENGTH50 = READ_LEFT_LENGTH50;

  private NgsParams buildParams(String leftReads, String rightReads) throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File leftok = FileUtils.createTempDir("left", "ngs", mDir);
    final File rightok = FileUtils.createTempDir("right", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(leftReads, leftok, null).close();
    ReaderTestUtils.getReaderDNA(rightReads, rightok, null).close();

    final NgsOutputParams op = NgsOutputParams.builder()
    .filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create()).create();
    return NgsParams.builder()
      .buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create())
      .buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create())
      .searchParams(SequenceParams.builder().directory(template).create())
      .maxFragmentLength(15).minFragmentLength(0).outputParams(op).create();
  }

  public void testConstructor() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(0);

    new SlidingWindowCollector(15, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));

    try {
      new SlidingWindowCollector(10, -1, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
      fail("Accepted bad minimum insert size");
    } catch (final IllegalArgumentException iae) {
      assertEquals("Minimum window size too small: -1", iae.getMessage());
    }

    try {
      new SlidingWindowCollector(0, 1, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
      fail("Accepted bad window size");
    } catch (final IllegalArgumentException iae) {
      assertEquals("Maximum window size too small: 0", iae.getMessage());
    }

    new SlidingWindowCollector(2, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));

    try {
      new SlidingWindowCollector(10, 0, MachineOrientation.ANY, null, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
      fail("Accepted null writer");
    } catch (final NullPointerException npe) {
      // expected
    }

    writer.close();
  }

  public void testNextTemplateId() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(4);

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);
    swc.nextTemplateId(2);
    swc.nextTemplateId(3);
    swc.nextTemplateId(4);

    writer.close(); // checks nextTemplateId call is passed through to writer
  }

  public void testMatch() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 6);

    final SlidingWindowCollector swc = new SlidingWindowCollector(11, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    //Two pairs, 0,4 and 4,13 (0,10 has too long a fragment length)
    swc.match(true, false, 0, 0);
    swc.match(false, true, 0, 4);
    swc.match(false, true, 0, 10);
    swc.match(true, false, 0, 13);

    swc.match(true, false, 1, 13);

    swc.nextTemplateId(2);

    writer.close();
    final MemoryPrintStream log = new MemoryPrintStream();
    Diagnostic.setLogStream(log.printStream());
    try {
      final Properties stats = swc.getStatistics();
      assertEquals("5", stats.getProperty("total_hits"));
      assertEquals("3", stats.getProperty("matings"));
      assertEquals("2", stats.getProperty("templates"));
      assertEquals(Integer.toString(11 + 64 + 64), stats.getProperty("window_size"));
      assertEquals("14", stats.getProperty("template_lengths_total"));
      final String results = log.toString();
      TestUtils.containsAll(results, "Reads with 0 potential right side candidates: 2 effective alignments: 0",
          "Reads with 1 potential right side candidates: 3 effective alignments: 6");
    } finally {
      log.close();
      Diagnostic.setLogStream();
    }
  }

  public void testMatchReadAtSameLocation() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 6);  //6 because 1,1 pairs both ways (4) + the 1,8 pair (2)

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    //two pairs, 1,1 and 1,8.
    swc.match(true, false, 0, 1);
    swc.match(false, true, 0, 1);
    swc.match(false, true, 0, 8);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testMatchNegative() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 4);

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 2, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    swc.match(true, false, 0, -9);
    swc.match(true, false, 0, -1);
    swc.match(false, true, 0, 0);
    swc.match(false, true, 0, 6);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testWindowPosition() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 4);
    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    //64 gets added to the windowlength
    assertEquals(0, swc.windowPosition(0));
    assertEquals(1, swc.windowPosition(1));
    assertEquals(73, swc.windowPosition(73));
    assertEquals(0, swc.windowPosition(138));
    assertEquals(1, swc.windowPosition(139));

    assertEquals(137, swc.windowPosition(-1));
    assertEquals(136, swc.windowPosition(-2));
    assertEquals(1, swc.windowPosition(-137));
    assertEquals(0, swc.windowPosition(-138));
    assertEquals(137, swc.windowPosition(-139));
    assertEquals(136, swc.windowPosition(-140));

    assertEquals(1, swc.windowPosition(-137 - 138));
    assertEquals(0, swc.windowPosition(-138 - 138));
    assertEquals(137, swc.windowPosition(-139 - 138));
    assertEquals(136, swc.windowPosition(-140 - 138));
  }

  public void testMatchReadAtSameLocationDifferentFrame() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 4);

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    swc.match(true, true, 0, 1);
    swc.match(true, false, 0, 1);
    swc.match(false, true, 0, 8);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testMatchWithLargeGap() throws Exception {
    final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 4);

    writer.setExpectedReadIds(new long[] {0, 0L, 0L, 0L});
    writer.setExpectedFirsts(new boolean[] {true, false, false, true});
    writer.setExpectedRC1s(new boolean[] {false, true, true, false});
    writer.setExpectedStart1s(new int[] {1, 8, 1011, 1017});

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    swc.match(true, false, 0, 1);
    swc.match(false, true, 0, 8);

    swc.match(false, true, 0, 1011);
    swc.match(true, false, 0, 1017);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testMatch3to3() throws Exception {
    final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 16);

    writer.setExpectedReadIds(new long[] {0, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L});
    writer.setExpectedFirsts(new boolean[] {true, true, true, true, true, true, true, true, false, false, false, false, false, false, false, false});
    writer.setExpectedRC1s(new boolean[] {false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true});
    //writer.setExpectedStart1s(new int[] {1, 1, 2, 2, 2, 3, 3, 3, 9, 9, 9, 10, 10, 10, 11, 11});
    writer.setExpectedStart1s(new int[]   {1, 1, 2, 2, 2, 3, 3, 3, 8, 8, 8, 9, 9, 9, 10, 10});

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    swc.match(true, false, 0, 1);
    swc.match(true, false, 0, 2);
    swc.match(true, false, 0, 3);

    swc.match(false, true, 0, 8);
    swc.match(false, true, 0, 9);
    swc.match(false, true, 0, 10);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testMatch3to3Longer() throws Exception {
    final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 0);
    final SlidingWindowCollector swc = new SlidingWindowCollector(400, 121, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH50, READ_RIGHT_LENGTH50));
    swc.nextTemplateId(1);
    swc.match(true, true, 0, 50);
    swc.match(false, false, 0, 120);
    swc.nextTemplateId(2);
    writer.close();
  }

  public void testMatch3to3Longer2() throws Exception {
    final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 2);
    final SlidingWindowCollector swc = new SlidingWindowCollector(400, 121, MachineOrientation.ANY, writer,  buildParams(READ_LEFT_LENGTH50, READ_RIGHT_LENGTH50));
    swc.nextTemplateId(1);
    swc.match(true, true, 0, 50);
    swc.match(false, false, 0, 220);
    swc.nextTemplateId(2);
    writer.close();
  }

  public void testMatch3to3Longer3() throws Exception {
    final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 0);
    final SlidingWindowCollector swc = new SlidingWindowCollector(400, 120, MachineOrientation.ANY, writer,  buildParams(READ_LEFT_LENGTH50, READ_RIGHT_LENGTH50));
    swc.nextTemplateId(1);
    swc.match(true, false, 0, 50);
    swc.match(false, true, 0, 401);
    swc.nextTemplateId(2);
    writer.close();
  }

  public void testMatchReadOutOfOrder() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 20);

    final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    swc.match(true, false, 1, 100);
    swc.match(false, true, 1, 108);

    swc.match(true, false, 1, 50);
    swc.match(true, false, 0, 50);
    swc.match(false, true, 1, 47);
    swc.match(false, true, 1, 46);
    swc.match(false, true, 0, 46);
    swc.match(false, true, 1, 44);
    swc.match(false, true, 0, 44);
    swc.match(false, true, 1, 45);
    swc.match(false, true, 0, 45);
    swc.match(false, true, 0, 47);
    swc.match(false, true, 0, 48);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testMatchMinSize() throws Exception {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 2);

    final SlidingWindowCollector swc = new SlidingWindowCollector(40, 20, MachineOrientation.ANY, writer, buildParams(READ_LEFT_LENGTH2, READ_RIGHT_LENGTH2));
    swc.nextTemplateId(1);

    swc.match(true, false, 0, 1);
    swc.match(false, true, 0, 8);
    swc.match(false, true, 0, 28);

    swc.nextTemplateId(2);

    writer.close();
  }

  public void testHitInfo() {
    final HitInfo hi = new HitInfo();
    hi.setValues(false, true, 1, 2);
    assertFalse(hi.first());
    assertTrue(hi.reverseComplement());
    assertEquals(1, hi.readId());
    assertEquals(2, hi.templateStart());
    assertEquals(-1, hi.score());
    assertFalse(hi.isPair(null));
    assertFalse(hi.isPair(hi));
    final HitInfo hi2 = new HitInfo();
    hi2.setValues(true, false, 1, 2);
    assertTrue(hi.isPair(hi2));
  }


  /**
   * Test a larger amount of data.
   * @throws Exception if there is a problem.
   */
  public void testMain() throws Exception {

    /** Discard all alignment output. */
    final PairedTempFileWriter samWriter = new PairedTempFileWriter() {
      @Override
      public boolean pairResultLeft(final MatedHitInfo matedHitInfo) {
        return true;
      }

      @Override
      public void pairResultRight(final MatedHitInfo matedHitInfo) {
      }

      @Override
      public void nextTemplateId(final long templateId) {
      }
      @Override
      public void close() {
      }
      @Override
      public void initialiseMated(OutputStream unmappedOut) {
      }
      @Override
      public void closeMated() {
      }

      @Override
      public void initialiseUnmated(OutputStream unmappedOut, MapQScoringReadBlocker a, MapQScoringReadBlocker b) {
      }
      @Override
      public boolean unmatedResult(int readId, boolean first, boolean rc, int start) {
        return true;
      }
      @Override
      public void closeUnmated() {
      }
    };

    final SlidingWindowCollector collector = new SlidingWindowCollector(80, 0, MachineOrientation.ANY, samWriter,  buildParams(READ_LEFT_LENGTH50, READ_RIGHT_LENGTH50));
    collector.nextTemplateId(1);
    int readId = 10;
    for (int i = 0; i < 10; i++) {
      collector.match(true, false, readId++, i);
      collector.match(true, false, readId++, i);
      collector.match(false, false, readId - 10, i + 1);
      collector.match(false, false, readId - 11, i + 1);
    }
    collector.nextTemplateId(Long.MAX_VALUE);

    final Properties stats = collector.getStatistics();
    //System.err.println(stats);
    assertEquals("40", stats.getProperty("total_hits"));
    assertEquals("2", stats.getProperty("templates"));
    assertEquals("11", stats.getProperty("template_lengths_total"));
    assertEquals(Integer.toString(80 + 64 + 64), stats.getProperty("window_size"));
    assertEquals("11", stats.getProperty("matings"));
  }

  static final String AS20 = "aaaaaaaaaaAAAAAAAAAA";
  static final String READ_LEFT_VARLENGTH = ">r1" + StringUtils.LS + AS50 + AS50 + StringUtils.LS;
  static final String READ_RIGHT_VARLENGTH = ">r1" + StringUtils.LS + AS20 + StringUtils.LS;


  public void testVarLength() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 2);

    final SlidingWindowCollector collector = new SlidingWindowCollector(200, 0, MachineOrientation.ANY, writer,  buildParams(READ_LEFT_VARLENGTH, READ_RIGHT_VARLENGTH));
    collector.nextTemplateId(1);
    collector.match(true, false, 0, 1);
    collector.match(false, true, 0, 181);
    collector.match(true, false, 0, 601);
    collector.match(false, true, 0, 782);
    collector.nextTemplateId(2);
    writer.close();

  }

  public void testDuplicateHitResolving() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 2);

    final SlidingWindowCollector collector = new SlidingWindowCollector(200, 0, MachineOrientation.ANY, writer,  buildParams(READ_LEFT_VARLENGTH, READ_RIGHT_VARLENGTH));
    collector.nextTemplateId(1);
    collector.match(true, false, 0, 1);
    collector.match(false, true, 2, 181);
    collector.match(false, true, 0, 181);
    collector.match(false, true, 1, 181);
    collector.match(false, true, 0, 181);
    collector.match(true, false, 0, 601);
    collector.match(false, true, 0, 782);
    collector.nextTemplateId(2);
    writer.close();
  }

  public void testDuplicateHitResolvingExtra() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 10);
    final SlidingWindowCollector collector = new SlidingWindowCollector(200, 0, MachineOrientation.ANY, writer,  buildParams(READ_LEFT_VARLENGTH, READ_RIGHT_VARLENGTH));
    collector.nextTemplateId(1);
    collector.match(true, false, 0, 1);
    collector.match(true, false, 0, 1);
    collector.match(true, false, 0, 2);
    collector.match(true, false, 0, 3);
    collector.match(true, false, 0, 3);
    collector.match(true, false, 0, 3);
    collector.match(true, false, 0, 3);
    collector.match(true, false, 0, 4);
    collector.match(true, false, 0, 5);
    collector.match(true, false, 0, 5);
    collector.match(true, false, 0, 3);
    collector.match(false, true, 0, 181);
    collector.nextTemplateId(2);
    writer.close();
  }

  public void testOppositeEndHash() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 2);

    final SlidingWindowCollector collector = new SlidingWindowCollector(200, 0, MachineOrientation.ANY, writer,  buildParams(READ_LEFT_VARLENGTH, READ_RIGHT_VARLENGTH));

    for (int i = 0; i < 10; i++) {
      assertEquals(i + (i % 2 == 0 ? 1 : -1), collector.readHash(i, true));
      assertEquals(i, collector.readHash(i, false));
    }
  }

  public void testNextPowerOfTwo() {
    for (int i = 0; i < 100000; i++) {
      int next = i;
      for (int j = 1; j <= 16; j *= 2) {
        next |= next >> j;
      }
      next++;
      assertEquals(next, SlidingWindowCollector.nextPowerOfTwo(i));
    }
  }

  private static final String READ_FARCICAL_TEST_1 = ">1\nacgt";

  public void testOrderProblemFarcical() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 2);
    final SlidingWindowCollector swc = new SlidingWindowCollector(30, 10, MachineOrientation.ANY, writer, buildParams(READ_FARCICAL_TEST_1, READ_FARCICAL_TEST_1));
    swc.nextTemplateId(0);
    swc.match(true, true, 0, 30);
    swc.match(false, false, 0, 114);
    swc.match(false, false, 0, 53);
    swc.nextTemplateId(Long.MAX_VALUE);
    writer.close();
  }

  private static final String READ_ORDER_TEST_1 = ">1\nacgtacgtacacgtacgtacacgtacgtacacgtacgtacacgtacgtacacgtacgtacacgtacgtacacgtacgtacacgtacgtacacgtacgtac";

  public void testOrderProblem() throws IOException {
    final PairedTempFileWriter writer = new TestPairedAlignmentWriter(2, 2);
    final SlidingWindowCollector swc = new SlidingWindowCollector(450, 400, MachineOrientation.ANY, writer, buildParams(READ_ORDER_TEST_1, READ_ORDER_TEST_1));
    swc.nextTemplateId(0);
    swc.match(true, true, 0, 72);
    //the next line should not mate with the first, but also should not cause flush preventing the next from mating
    swc.match(false, false, 0, 430);
    swc.match(false, false, 0, 393);
    swc.nextTemplateId(Long.MAX_VALUE);
    writer.close();
  }


  private void addRandomRead(StringBuilder readsLeft, int i, char[] chars, Random rand) {
    readsLeft.append(">r");
    readsLeft.append(i);
    readsLeft.append(StringUtils.LS);
    for (int j = 0; j < 36; j++) {
      readsLeft.append(chars[rand.nextInt(4)]);
    }
    readsLeft.append(StringUtils.LS);
  }

  //Many reads;
  public void testMaxHitsInPosition() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 2);
      final StringBuilder readsLeft = new StringBuilder();
      final StringBuilder readsRight = new StringBuilder();
      final char[] chars = {'a', 'c', 'g', 't'};
      final Random rand = new Random(33);
      final int maxHits = GlobalFlags.getIntegerValue(GlobalFlags.SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG);
      for (int i = 0; i < maxHits + 2; i++) {
        addRandomRead(readsLeft, i, chars, rand);
        addRandomRead(readsRight, i, chars, rand);
      }

      final SlidingWindowCollector swc = new SlidingWindowCollector(45, 0, MachineOrientation.ANY, writer, buildParams(readsLeft.toString(), readsRight.toString()));
      swc.nextTemplateId(0);

      // add 1000 left arms at a position which should be ignored
      for (int i = 0; i < maxHits + 1; i++) {
        swc.match(true, false, i, 100);
      }
      swc.match(true, false, 0, 99);

      // add 1000 right arms at a position which should be ignored
      for (int i = 0; i < maxHits + 1; i++) {
        swc.match(false, false, i, 108);
      }
      swc.match(false, false, 0, 107);

      swc.nextTemplateId(1);

      writer.close();
      TestUtils.containsAll(mps.toString(),
            "Max hits per position exceeded at template: 0 templateStart: 100",
            "Max hits per position exceeded at template: 0 templateStart: 108");

    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testMaxMatedHitsInPosition() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final int[] count = new int[1];
    try {
      final TestPairedAlignmentWriter writer = new TestPairedAlignmentWriter(2, 2) {
        @Override
        public void pairResultRight(final MatedHitInfo matedHitInfo) {
          count[0]++;
        }
      };
      final StringBuilder readsLeft = new StringBuilder();
      final StringBuilder readsRight = new StringBuilder();

      final char[] chars = {'a', 'c', 'g', 't'};
      final Random rand = new Random(33);
      final int maxHits = GlobalFlags.getIntegerValue(GlobalFlags.SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG);
      for (int i = 0; i < maxHits + 2; i++) {
        addRandomRead(readsLeft, i, chars, rand);
        addRandomRead(readsRight, i, chars, rand);
      }

      final SlidingWindowCollector swc = new SlidingWindowCollector(10, 0, MachineOrientation.ANY, writer, buildParams(readsLeft.toString(), readsRight.toString()));
      swc.nextTemplateId(0);

      for (int i = 0; i < maxHits + 1; i++) {
        final HitInfo hit = new HitInfo();
        final HitInfo mate = new HitInfo();
        mate.mTemplateStart = 108;

        swc.checkPair(hit, mate);

      }
      assertEquals(maxHits, writer.mPairedLeft); //we know we've called checkPair max hits + 1 times, but it should have only called pairResultLeft max hits times.
      TestUtils.containsAll(mps.toString(), "Max hits per position exceeded at template: 0 templateStart: unknown");
      swc.mCurrentReferencePosition = 100;
      swc.flushToPosition(110);

      assertEquals(1000, count[0]);
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
