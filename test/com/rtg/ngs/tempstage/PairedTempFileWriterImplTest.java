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
package com.rtg.ngs.tempstage;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.ngs.ReadStatusTracker;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.pairedend.MatedHitInfo;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests corresponding class.
 */
public class PairedTempFileWriterImplTest extends TestCase {

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

  /**
   * a useless status ID listener
   */
  public static class UselessStatusIdListener implements ReadStatusListener {
    @Override
    public void addStatus(int readId, int status) {
    }

    @Override
    public String toString() {
      return "UselessStatusIdListener";
    }
  }

  public PairedTempFileWriterImpl createPairedWriter(final NgsParams param, final File outputFile, boolean writeHeader) throws IOException {
    final PairedTempFileWriterImpl impl = new PairedTempFileWriterImpl(param, new UselessStatusIdListener(), SharedResources.generateSharedResources(param));
    impl.initialiseMated(FileUtils.createOutputStream(outputFile, false, false));
    return impl;
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg" + StringUtils.LS;
  static final String TEMP_LEFT = "acacactgcaagcaagagggcctccc"; //1 deletion vs TEMPLATE
  static final String TEMP_RIGHT = "cccctttggcccccgaccagtgtgggctga"; //perfect match, reverse compliment
  static final String READ_LEFT = ">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS;
  static final String READ_RIGHT = ">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS;

  MatedHitInfo populate(final MatedHitInfo mh, final int readId, final boolean first, final boolean rc1, final int start1, final boolean rc2, final int start2) {
    mh.setValues(readId, !first, rc1, start1, rc2, start2);
    return mh;
  }

  public void testPairing() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);
    final File out = File.createTempFile("sam", "out", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

    CommandLine.setCommandArgs("wibble", "-h", "dream-turnip");
    try {
      final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create()).create();
      final NgsParams param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left)
        .useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create())
        .maxFragmentLength(1000).minFragmentLength(0).outputParams(op)
        .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0)
        .create();

      final PairedTempFileWriterImpl w = createPairedWriter(param, out, false);
      try {
        final MatedHitInfo mh = new MatedHitInfo();
        try {
          w.pairResultLeft(mh);
          fail();
        } catch (final RuntimeException re) {
          assertEquals("Call nextTemplateId first", re.getMessage());
        }
        w.nextTemplateId(0);
        if (w.pairResultLeft(populate(mh, 0, true, false, 0, true, 31))) {
          w.pairResultRight(mh);
        }
        w.nextTemplateId(Long.MAX_VALUE);
        try {
          w.pairResultLeft(populate(mh, 0, false, true, 0, true, 0));
          fail();
        } catch (final RuntimeException re) {
          assertEquals("Call nextTemplateId first", re.getMessage());
        }
      } finally {
        w.close();
        w.closeMated();
      }
      assertNotNull(w.getBlocker());

      final TempRecordReaderNio dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
      try {
        BinaryTempFileRecord rec = dis.readRecord();
        assertNotNull(rec);
        checkRecord(rec, 0, 99, 0, 1, "12=1D14=", 37, 66, 2);
        rec = dis.readRecord();
        assertNotNull(rec);
        checkRecord(rec, 0, 147, 0, 37, "30=", 1, -66, 0);
        assertNull(dis.readRecord());
      } finally {
        dis.close();
      }
    } finally {
      CommandLine.clearCommandArgs();
    }
  }

  static final String TEMPLATE2 = TEMPLATE + ">x" + StringUtils.LS + "ttttttttttttttttttttttttttttttttttttttttttttttt" + StringUtils.LS;

  public void testDualTemplateReverseOrder() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);
    final File out = File.createTempFile("sam", "out", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE2, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();

    final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create()).create();
    final NgsParams param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0).outputParams(op).create();

    final PairedTempFileWriterImpl w = createPairedWriter(param, out, true);
    final MatedHitInfo mh = new MatedHitInfo();
    try {
      w.nextTemplateId(1);
      w.pairResultLeft(populate(mh, 0, false, true, 26, false, 0));
      try {
        w.nextTemplateId(0);
        fail();
      } catch (final IllegalArgumentException e) {
        // ok
      }
    } finally {
      w.close();
    }
  }

  static final String TEMP_LEFT_MED = "tgcaaccaagagggcctccc";  //1 substitutions
  static final String TEMP_RIGHT_MED = "tttggcccccgaccagtctgggctga"; //1 substitutions
  static final String TEMP_LEFT_BAD = "tgcaaccatgaggccctccc";  //3 substitutions
  static final String TEMP_RIGHT_BAD = "tttgacctccgaccagtctgggctga"; //3 substitutions
  static final String READ_LEFT_BAD = ">t" + StringUtils.LS + TEMP_LEFT_BAD + StringUtils.LS;  //3 substitutions
  static final String READ_RIGHT_BAD = ">t" + StringUtils.LS + TEMP_RIGHT_BAD + StringUtils.LS; //3 substitutions
  static final String TEMPLATE_3 = ">t" + StringUtils.LS + "tttaccccccccccccc" + StringUtils.LS;
  static final String TEMP_LEFT_3 = "ttt";
  static final String TEMP_RIGHT_3 = "ccc";
  static final String READ_LEFT_3 = ">r" + StringUtils.LS + TEMP_LEFT_3 + StringUtils.LS;
  static final String READ_RIGHT_3 = ">r" + StringUtils.LS + TEMP_RIGHT_3 + StringUtils.LS;

  public void testMaxAlignmentScore() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File leftok = FileUtils.createTempDir("left", "ngs", mDir);
    final File rightok = FileUtils.createTempDir("right", "ngs", mDir);
    final File leftbad = FileUtils.createTempDir("leftbad", "ngs", mDir);
    final File rightbad = FileUtils.createTempDir("rightbad", "ngs", mDir);
    final File templatetrio = FileUtils.createTempDir("templatetrio", "ngs", mDir);
    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, leftok, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, rightok, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT_BAD, leftbad, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT_BAD, rightbad, null).close();
    ReaderTestUtils.getReaderDNA(">t" + StringUtils.LS + TEMP_LEFT_MED + TEMP_RIGHT_MED + TEMP_LEFT_BAD + TEMP_RIGHT_BAD + TEMP_LEFT + TEMP_RIGHT + "nnnnnnn" + StringUtils.LS, templatetrio, null).close();

    NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().matedMaxMismatches(IntegerOrPercentage.valueOf(2)).outputFilter(OutputFilter.PAIRED_END).create()).create();
    NgsParams param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(rightbad).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0).outputParams(op).create();
    final MockStatusTracker sril = new MockStatusTracker();
    SharedResources sr = SharedResources.generateSharedResources(param);
    PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, sril, sr);
    final MatedHitInfo mh = new MatedHitInfo();
    final MemoryPrintStream stream = new MemoryPrintStream();
    try {
      w.initialiseMated(stream.outputStream());
      w.nextTemplateId(0);
      if (w.pairResultLeft(populate(mh, 0, true, false, 0, true, 26))) {
        w.pairResultRight(mh);
      }
    } finally {
      w.close();
    }
    final byte[] actual = IOUtils.readData(new ByteArrayInputStream(stream.toByteArray()));
    assertTrue(Arrays.equals(new byte[]{(byte) 0xff, (byte) 0xff, (byte) 0xff, (byte) 0xff}, actual));
    assertEquals(",0:MATED ", sril.mStatusString);
    sril.mStatusString = "";

    param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(leftbad).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0).outputParams(op).create();
    sr = SharedResources.generateSharedResources(param);
    w = new PairedTempFileWriterImpl(param, sril, sr);
    try {
      w.nextTemplateId(0);
      if (w.pairResultLeft(populate(mh, 0, true, false, 0, true, 26))) {
        w.pairResultRight(mh);
      }
    } finally {
      w.close();
    }
    assertEquals(",0:MATED ", sril.mStatusString);
    sril.mStatusString = "";

    param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(leftbad).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(rightbad).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0).outputParams(op).create();
    sr = SharedResources.generateSharedResources(param);
    w = new PairedTempFileWriterImpl(param, sril, sr);
    try {
      w.nextTemplateId(0);
      if (w.pairResultLeft(populate(mh, 0, true, false, 0, true, 26))) {
        w.pairResultRight(mh);
      }
    } finally {
      w.close();
    }
    assertEquals(",0:MATED ", sril.mStatusString);
    sril.mStatusString = "";

    op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().matedMaxMismatches(IntegerOrPercentage.valueOf(2)).outputFilter(OutputFilter.PAIRED_END).create()).create();
    param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create())
      .maxFragmentLength(1000).minFragmentLength(0).outputParams(op)
      .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0)
      .create();
    sr = SharedResources.generateSharedResources(param);
    w = new PairedTempFileWriterImpl(param, sril, sr);
    final File out = File.createTempFile("sam", "out", mDir);
    final FileOutputStream toBeSure = new FileOutputStream(out);
    try {
      try {
        w.initialiseMated(toBeSure);
        w.nextTemplateId(0);
        if (w.pairResultLeft(populate(mh, 0, true, false, 0, true, 31))) {
          w.pairResultRight(mh);
        }
      } finally {
        w.close();
      }
    } finally {
      toBeSure.close();
    }
    TempRecordReaderNio dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 99, 0, 1, "12=1D14=", 37, 66, 2);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 147, 0, 37, "30=", 1, -66, 0);
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }

    assertEquals(",0:MATED ,0:MATED_ALIGN_SCORE ,0:MATED ,0:MATED_ALIGN_SCORE ", sril.mStatusString);
    sril.mStatusString = "";

    op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().matedMaxMismatches(IntegerOrPercentage.valueOf(20)).outputFilter(OutputFilter.PAIRED_END).create()).create();
    param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(templatetrio).loadNames(true).useMemReader(true).create())
      .maxFragmentLength(1000).minFragmentLength(0).outputParams(op).alignerBandWidthFactor(new MaxShiftFactor(0.55))
      .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(1)
      .create();
    w = createPairedWriter(param, out, true);
    try {
      w.nextTemplateId(0);
      if (w.pairResultLeft(populate(mh, 0, true, false, 0, true, 26))) {
        w.pairResultRight(mh);
      }
      if (w.pairResultLeft(populate(mh, 0, true, false, 46, true, 72))) {
        w.pairResultRight(mh);
      }
      if (w.pairResultLeft(populate(mh, 0, true, false, 56, true, 82))) {
        w.pairResultRight(mh);
      }
      if (w.pairResultLeft(populate(mh, 0, true, false, 93, true, 120))) {
        w.pairResultRight(mh);
      }
    } finally {
      w.close();
    }

    dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 99, 0, 1, "6S5=1X14=", 28, 63, 7);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 147, 0, 28, "1X1=1X1=1X3=2X4=1X1=1X1=2X2=1X1=4X1=1X", 1, -63, 15);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 99, 0, 42, "1X1=2X1=1I5=1X2=1X4=1X6=", 74, 62, 8);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 147, 0, 74, "2=1X1=1X3=2X4=1X1=1X1=4X2=2X1=3X", 42, -62, 15);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 99, 0, 93, "26=", 115, 51, 0);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 147, 0, 115, "2=2X3=1I1=2X3=1X1=3X1=1X3=2X1=1X2=", 93, -51, 14);
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }
  }

  public void testReadBlocking() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);
    final File out = File.createTempFile("sam", "out", mDir);

    ReaderTestUtils.getReaderDNA(READ_LEFT_3, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT_3, right, null).close();

    final String templ = TEMPLATE_3;

    ReaderTestUtils.getReaderDNA(templ, template, null).close();

    final NgsParams param = getCommonTestParams(left, right, template, IntegerOrPercentage.valueOf(5), IntegerOrPercentage.valueOf(5));

    final PairedTempFileWriterImpl w = createPairedWriter(param, out, true);
    final MatedHitInfo mh = new MatedHitInfo();
    try {
      try {
        w.pairResultLeft(populate(mh, 0, false, true, 0, true, 0));
        fail();
      } catch (final RuntimeException re) {
        assertEquals("Call nextTemplateId first", re.getMessage());
      }
      pairResults3(w);
      w.nextTemplateId(Long.MAX_VALUE);
    } finally {
      w.close();
    }

    final TempRecordReaderNio dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      BinaryTempFileRecord rec1;
      for (int i = 6; i <= 15; i++) {
        rec1 = dis.readRecord();
        assertNotNull(rec1);
        checkRecord(rec1, 0, 67, 0, 1, "3=", i, i + 2, 0);
      }
      for (int i = 6; i <= 15; i++) {
        rec1 = dis.readRecord();
        assertNotNull(rec1);
        checkRecord(rec1, 0, 131, 0, i, "3=", 1, -(i + 2), 0);
      }
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }
  }

  private void pairResults3(final PairedTempFileWriterImpl w) throws IOException {
    w.nextTemplateId(0);
    final MatedHitInfo[] mhs = new MatedHitInfo[11];
    for (int i = 0; i < mhs.length; i++) {
      final MatedHitInfo mh = new MatedHitInfo();
      if (w.pairResultLeft(populate(mh, 0, true, false, 0, false, i + 5))) {
        mhs[i] = mh;
      }
    } // Pair on final iteration should be blocked

    for (final MatedHitInfo mh : mhs) {
      if (mh != null) {
        w.pairResultRight(mh);
      }
    } // Pair on final iteration should be blocked
  }

  public void testCheckPairScores() throws Exception {
    final File out = File.createTempFile("sam", "out", mDir);
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(READ_LEFT, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, right, null).close();
    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();

    CommandLine.setCommandArgs("wibble", "-h", "dream-turnip");
    try {
      final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create()).create();
      final NgsParams param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create())
        .buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create())
        .searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0)
        .outputParams(op)
        .substitutionPenalty(3)
        .create();

      try (PairedTempFileWriterImpl w = createPairedWriter(param, out, true)) {
        w.nextTemplateId(0);
        assertTrue(w.checkPairScores(0, false, 9, 2));
        assertFalse(w.checkPairScores(0, false, 10, 2));
        assertTrue(w.checkPairScores(0, false, 3, 6));
        assertFalse(w.checkPairScores(0, false, 3, 7));

        assertFalse(w.checkPairScores(0, false, 10, 7));
        final MemoryPrintStream mps = new MemoryPrintStream();
        try {
          Diagnostic.setLogStream(mps.printStream());
          w.closeMated();
          assertTrue(mps.toString().contains("(mismatch threshold of 10%): 2 passed, 5 total"));
        } finally {
          Diagnostic.setLogStream();
        }
      }
    } finally {
      CommandLine.clearCommandArgs();
    }
  }

  private static NgsParams getCommonTestParams(File left, File right, File template, IntegerOrPercentage maxMatedScore, IntegerOrPercentage maxUnmatedScore) {
    final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder()
      .matedMaxMismatches(maxMatedScore)
      .unmatedMaxMismatches(maxUnmatedScore)
      .outputFilter(OutputFilter.PAIRED_END).maxTopResults(10).create()).create();
    return NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create())
      .maxFragmentLength(1000).minFragmentLength(0).outputParams(op)
      .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(1)
      .create();
  }

  public void testOtherStuff() throws Exception {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);
    final File out = File.createTempFile("sam", "out", mDir);

    ReaderTestUtils.getReaderDNA(READ_LEFT_3, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT_3, right, null).close();

    final String templ = TEMPLATE_3;

    ReaderTestUtils.getReaderDNA(templ, template, null).close();
    final NgsParams param = getCommonTestParams(left, right, template, IntegerOrPercentage.valueOf(5), IntegerOrPercentage.valueOf(5));


    final MemoryPrintStream baos = new MemoryPrintStream();
    final ByteArrayOutputStream baosunmated = new ByteArrayOutputStream();
    Diagnostic.setLogStream(baos.printStream());
    try {
      final SharedResources sr = SharedResources.generateSharedResources(param);
      final ReadStatusListener sril = new UselessStatusIdListener();
      final PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, sril, sr);
      final FileOutputStream toBeSure = new FileOutputStream(out);
      try {
        w.initialiseMated(toBeSure);
        try {
          pairResults3(w);

          try {
            w.unmatedResult(0, true, true, 0);
            fail();
          } catch (final NullPointerException npe) {
            //expected - must call initialiseUnmated first!
          }

          final MapQScoringReadBlocker srbl = new MapQScoringReadBlocker(1, 1);
          final MapQScoringReadBlocker srbr = new MapQScoringReadBlocker(1, 1);
          w.initialiseUnmated(baosunmated, srbl, srbr);
          try {
            w.nextTemplateId(0);
            w.unmatedResult(0, true, true, 0);
            assertEquals(1, srbl.getCount1(0));
            assertEquals(0, srbr.getCount1(0));
            w.nextTemplateId(Long.MAX_VALUE);
          } finally {
            w.closeUnmated();
            srbl.close();
            srbr.close();
          }
        } finally {
          w.close();
          assertNull(w.mTemplateReader);
          assertNull(w.mSecondReader);
          assertNull(w.mFirstReader);
        }
      } finally {
        toBeSure.close();
      }


      final String diaglogstring = baos.toString();
      assertTrue(diaglogstring.contains("Statistics of blocked pairings"));
      assertTrue(diaglogstring, diaglogstring.contains("Duplicates detected during SAM writing: 2"));
      assertTrue(diaglogstring, diaglogstring.contains("Total reads 1"));
      assertTrue(diaglogstring, diaglogstring.contains("1 reads had count 10"));
      assertFalse(diaglogstring, diaglogstring.contains("0\t89\tt"));
      assertFalse(diaglogstring, diaglogstring.contains("0\t77\t*"));
      assertFalse(diaglogstring, diaglogstring.contains("0\t67\tt"));
      assertFalse(diaglogstring, diaglogstring.contains("0\t131\tt"));


      TempRecordReaderNio dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
      try {
        BinaryTempFileRecord rec1;
        for (int i = 6; i <= 15; i++) {
          rec1 = dis.readRecord();
          assertNotNull(rec1);
          checkRecord(rec1, 0, 67, 0, 1, "3=", i, i + 2, 0);
        }
        for (int i = 6; i <= 15; i++) {
          rec1 = dis.readRecord();
          assertNotNull(rec1);
          checkRecord(rec1, 0, 131, 0, i, "3=", 1, -(i + 2), 0);
        }
        assertNull(dis.readRecord());
      } finally {
        dis.close();
      }

      dis = new TempRecordReaderNio(new ByteArrayInputStream(baosunmated.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, false));
      try {
        final BinaryTempFileRecord rec1 = dis.readRecord();
        assertNotNull(rec1);
        checkRecord(rec1, 0, 81, 0, 2, "2X1=", 0, 0, 2);
        assertNull(dis.readRecord());
      } finally {
        dis.close();
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testUselessStatusIdListener() {
    final UselessStatusIdListener listener = new UselessStatusIdListener();
    assertEquals(listener.toString(), "UselessStatusIdListener");

  }

  NgsParams unmatedParams() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File leftok = FileUtils.createTempDir("left", "ngs", mDir);
    final File rightok = FileUtils.createTempDir("right", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, leftok, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, rightok, null).close();


    final NgsOutputParams op = NgsOutputParams.builder()
      .filterParams(NgsFilterParams.builder().matedMaxMismatches(IntegerOrPercentage.valueOf(2))
        .outputFilter(OutputFilter.PAIRED_END).create())
      .create();
    return NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create()).buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create()).searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create())
      .maxFragmentLength(1000).minFragmentLength(0).outputParams(op)
      .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0)
      .create();
  }

  public void testUnmatedUnfiltered() throws IOException {

    final NgsParams param = unmatedParams();

    final SharedResources sr = SharedResources.generateSharedResources(param);
    final MockStatusTracker sril = new MockStatusTracker();

    try (PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, sril, sr)) {
      final ByteArrayOutputStream baosunmated = new ByteArrayOutputStream();
      final MapQScoringReadBlocker srbl = new MapQScoringReadBlocker(1, 1);
      final MapQScoringReadBlocker srbr = new MapQScoringReadBlocker(1, 1);
      w.initialiseUnmated(baosunmated, srbl, srbr);

      w.nextTemplateId(0);
      final BinaryTempFileRecord sam = w.unmatedResultUnfiltered(0, true, false, 0);

      assertEquals(0, sam.getMatePosition());
      assertEquals(false, sam.isReverseStrand());

      assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_ALIGN_SCORE_FIRST ,0:UNMATED_FIRST ", sril.mStatusString);
      sril.mStatusString = "";

      final BinaryTempFileRecord samRight = w.unmatedResultUnfiltered(0, false, true, 38);
      assertNotNull(samRight);
      assertEquals(0, samRight.getMatePosition());
      assertTrue(samRight.isReverseStrand());
      assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_SECOND ,0:UNMATED_ALIGN_SCORE_SECOND ,0:UNMATED_SECOND ", sril.mStatusString);
    }
  }

  public void testUnmated() throws IOException {

    final NgsParams param = unmatedParams();

    final SharedResources sr = SharedResources.generateSharedResources(param);
    final MockStatusTracker sril = new MockStatusTracker();

    final ByteArrayOutputStream baosunmated = new ByteArrayOutputStream();
    try (PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, sril, sr)) {
      final MapQScoringReadBlocker srbl = new MapQScoringReadBlocker(1, 1);
      final MapQScoringReadBlocker srbr = new MapQScoringReadBlocker(1, 1);
      w.initialiseUnmated(baosunmated, srbl, srbr);

      w.nextTemplateId(0);
      boolean result = w.unmatedResult(0, true, false, 0);

      assertTrue(result);

      assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_ALIGN_SCORE_FIRST ", sril.mStatusString);
      sril.mStatusString = "";

      result = w.unmatedResult(0, false, true, 40);
      assertTrue(result);

      assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_SECOND ,0:UNMATED_ALIGN_SCORE_SECOND ", sril.mStatusString);
      sril.mStatusString = "";
    }
    final TempRecordReaderNio dis = new TempRecordReaderNio(new ByteArrayInputStream(baosunmated.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 65, 0, 1, "12=1D14=", 0, 0, 2);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 145, 0, 37, "30=", 0, 0, 0);
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }
  }

  public void testUnmatedRegion() throws IOException {

    final NgsParams param = unmatedParams();

    final SharedResources sr = SharedResources.generateSharedResources(param);
    final MockStatusTracker sril = new MockStatusTracker();
    final PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, sril, sr);

    final ByteArrayOutputStream baosunmated = new ByteArrayOutputStream();
    w.setClipRegion(new HashingRegion(0, 0, 0, 20, 0, 20));
    try {
      final MapQScoringReadBlocker srbl = new MapQScoringReadBlocker(1, 1);
      final MapQScoringReadBlocker srbr = new MapQScoringReadBlocker(1, 1);
      w.initialiseUnmated(baosunmated, srbl, srbr);

      w.nextTemplateId(0);
      boolean result = w.unmatedResult(0, true, false, 0);

      assertTrue(result);

      assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_FIRST ,0:UNMATED_ALIGN_SCORE_FIRST ", sril.mStatusString);
      sril.mStatusString = "";

      result = w.unmatedResult(0, false, true, 40);
      assertTrue(result);

      assertEquals(",0:UNMATED_COMPUTE_ALIGNMENT_SECOND ", sril.mStatusString);
      sril.mStatusString = "";


    } finally {
      w.close();
    }
    final TempRecordReaderNio dis = new TempRecordReaderNio(new ByteArrayInputStream(baosunmated.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      final BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 65, 0, 1, "12=1D14=", 0, 0, 2);
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }
  }


  public void testUnmatedUnfilteredRegion() throws IOException {
    final NgsParams param = unmatedParams();

    final SharedResources sr = SharedResources.generateSharedResources(param);
    final MockStatusTracker sril = new MockStatusTracker();
    final PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, sril, sr);
    w.setClipRegion(new HashingRegion(0, 0, 0, 20, 0, 20));

    try {
      final ByteArrayOutputStream baosunmated = new ByteArrayOutputStream();
      final MapQScoringReadBlocker srbl = new MapQScoringReadBlocker(1, 1);
      final MapQScoringReadBlocker srbr = new MapQScoringReadBlocker(1, 1);
      w.initialiseUnmated(baosunmated, srbl, srbr);

      w.nextTemplateId(0);
      final BinaryTempFileRecord sam = w.unmatedResultUnfiltered(0, true, false, 0);

      assertEquals(0, sam.getMatePosition());
      assertEquals(false, sam.isReverseStrand());

      final BinaryTempFileRecord samRight = w.unmatedResultUnfiltered(0, false, true, 40);
      assertNull(samRight);

      w.setClipRegion(HashingRegion.NONE);
      final BinaryTempFileRecord samRight2 = w.unmatedResultUnfiltered(0, false, true, 40);
      assertNotNull(samRight2);

    } finally {
      w.close();
    }
  }

  static class MockStatusTracker implements ReadStatusListener {

    protected String mStatusString = "";

    @Override
    public void addStatus(int readId, int status) {
      mStatusString = mStatusString + "," + readId + ":" + ReadStatusTracker.statusToString(status);
    }
  }

  public void testPairRegion() throws Exception {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);
    final File out = File.createTempFile("sam", "out", mDir);

    ReaderTestUtils.getReaderDNA(READ_LEFT_3, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT_3, right, null).close();

    final String templ = TEMPLATE_3;

    ReaderTestUtils.getReaderDNA(templ, template, null).close();
    final NgsParams param = getCommonTestParams(left, right, template, IntegerOrPercentage.valueOf(5), IntegerOrPercentage.valueOf(5));

    try {
      final PairedTempFileWriterImpl w = createPairedWriter(param, out, true);
      try {
        w.setClipRegion(new HashingRegion(0, 6, 0, 10, 7, 11));
        pairResults3(w);

      } finally {
        w.close();
        assertNull(w.mTemplateReader);
        assertNull(w.mSecondReader);
        assertNull(w.mFirstReader);
      }
    } finally {
      Diagnostic.setLogStream();
    }
    final TempRecordReaderNio dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 131, 0, 7, "3=", 1, null, null);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 131, 0, 8, "3=", 1, null, null);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 131, 0, 9, "3=", 1, null, null);
      rec = dis.readRecord();
      assertNotNull(rec);
      checkRecord(rec, 0, 131, 0, 10, "3=", 1, null, null);
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }


  }

  void checkRecord(BinaryTempFileRecord rec, int readId, int samFlags, int refId, int startPos, String cigar, int mateStartPos, Integer tlen, Integer as) {
    assertEquals(readId, rec.getReadId());
    assertEquals(samFlags, rec.getSamFlags() & 0xff);
    assertEquals(refId, rec.getReferenceId());
    assertEquals(startPos, rec.getStartPosition());
    assertEquals(cigar, new String(rec.getCigarString()));
    assertEquals(mateStartPos, rec.getMatePosition());
    if (tlen != null) {
      assertEquals(tlen.intValue(), rec.getTemplateLength());
    }
    if (as != null) {
      assertEquals(as.intValue(), rec.getAlignmentScore());
    }
  }

  public void testTemplateOffset() throws Exception {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);
    final File out = File.createTempFile("sam", "out", mDir);
    final File unmated = File.createTempFile("sam", "unmated", mDir);

    ReaderTestUtils.getReaderDNA(READ_LEFT_3, left, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT_3, right, null).close();

    final StringBuilder templateBuilder = new StringBuilder();
    templateBuilder.append(">t0").append(StringUtils.LS);
    for (int i = 0; i < 1500; i++) {
      templateBuilder.append("a");
    }
    templateBuilder.append("tttaccccccccccccc");
    for (int i = 0; i < 1500; i++) {
      templateBuilder.append("a");
    }
    ReaderTestUtils.getReaderDNA(templateBuilder.toString(), template, null).close();
    final NgsParams param = getCommonTestParams(left, right, template, IntegerOrPercentage.valueOf(5), IntegerOrPercentage.valueOf(5));
    final MapQScoringReadBlocker srbl = new MapQScoringReadBlocker(1, 1);
    final MapQScoringReadBlocker srbr = new MapQScoringReadBlocker(1, 1);
    try {
      final PairedTempFileWriterImpl w = createPairedWriter(param, out, true);
      try {
        w.initialiseUnmated(FileUtils.createOutputStream(unmated, false), srbl, srbr);
        w.setClipRegion(new HashingRegion(0, 1506, 0, 1550, 1506, 1550));
        w.nextTemplateId(0);
        w.unmatedResultUnfiltered(0, false, true, 1505);
        w.unmatedResultUnfiltered(0, false, true, 1506);
        final MatedHitInfo first = new MatedHitInfo();
        final MatedHitInfo second = new MatedHitInfo();
        w.pairResultLeft(populate(first, 0, true, false, 0, false, 1505));
        w.pairResultLeft(populate(second, 0, true, false, 0, false, 1506));

        w.pairResultRight(first);
        w.pairResultRight(second);

      } finally {
        w.close();
        assertNull(w.mTemplateReader);
        assertNull(w.mSecondReader);
        assertNull(w.mFirstReader);
      }
    } finally {
      Diagnostic.setLogStream();
    }
    TempRecordReaderNio dis = new TempRecordReaderNio(new FileInputStream(out), new TempRecordReader.RecordFactory(true, false, false, false));
    try {
      final BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      assertEquals(0, rec.getReadId());
      assertEquals(131, rec.getSamFlags() & 0xff);
      assertEquals(0, rec.getReferenceId());
      assertEquals(1507, rec.getStartPosition());
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }

    dis = new TempRecordReaderNio(new FileInputStream(unmated), new TempRecordReader.RecordFactory(true, false, false, true));
    try {
      final BinaryTempFileRecord rec = dis.readRecord();
      assertNotNull(rec);
      assertEquals(0, rec.getReadId());
      assertEquals(145, rec.getSamFlags() & 0xff);
      assertEquals(0, rec.getReferenceId());
      assertEquals(1507, rec.getStartPosition());
      assertNull(dis.readRecord());
    } finally {
      dis.close();
    }
  }

  public void testAlignmentArm() throws Exception {
    final File out = File.createTempFile("sam", "out", mDir);
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File left = FileUtils.createTempDir("left", "ngs", mDir);
    final File right = FileUtils.createTempDir("right", "ngs", mDir);

    ReaderTestUtils.getReaderDNAFastqCG("@s1\nctggttaaaatatgaagtgaccaccatgcttgaga\n+\nabcdefghijklmnopqrstuvwxyz123456789\n", left, PrereadArm.LEFT).close();
    ReaderTestUtils.getReaderDNAFastqCG("@s1\ntctcaagcatggtggtcacttcatattttaaccag\n+\nabcdefghijklmnopqrstuvwxyz123456789\n", right, PrereadArm.RIGHT).close();
    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();

    CommandLine.setCommandArgs("wibble", "-h", "dream-turnip");
    try {
      final NgsOutputParams op = NgsOutputParams.builder().filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create()).create();
      final NgsParams param = NgsParams.builder().buildFirstParams(SequenceParams.builder().directory(left).useMemReader(true).create())
        .buildSecondParams(SequenceParams.builder().directory(right).useMemReader(true).create())
        .searchParams(SequenceParams.builder().directory(template).loadNames(true).useMemReader(true).create()).maxFragmentLength(1000).minFragmentLength(0)
        .substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1).unknownsPenalty(0)
        .outputParams(op).create();

      final PairedTempFileWriterImpl w = new PairedTempFileWriterImpl(param, new UselessStatusIdListener(), SharedResources.generateSharedResources(param)) {
        @Override
        protected final int[] calculateEditDistance(byte[] read, int length, int start, boolean rc, IntegerOrPercentage maxMismatches, boolean first, int readId) {
          assertFalse(first);
          return super.calculateEditDistance(read, length, start, rc, maxMismatches, first, readId);
        }
      };
      w.initialiseMated(FileUtils.createOutputStream(out, false, false));
      try {
        w.nextTemplateId(0);
        final MatedHitInfo mhi = new MatedHitInfo();
        mhi.setValues(0, true, false, 0, false, 0);
        w.pairResultLeft(mhi);
      } finally {
        w.close();
        w.closeMated();
      }
    } finally {
      CommandLine.clearCommandArgs();
    }
  }
}
