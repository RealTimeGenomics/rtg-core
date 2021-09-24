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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.zip.GZIPOutputStream;

import com.rtg.bed.BedUtils;
import com.rtg.calibrate.CalibratorTest;
import com.rtg.mode.SequenceType;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.TempRecordWriter;
import com.rtg.ngs.tempstage.TempRecordWriterNio;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import junit.framework.TestCase;

/**
 */
public class DummySamResultsFilterTest extends TestCase {

  private final class DummySamResultsFilter extends AbstractSamResultsFilter {
    private final MapQScoringReadBlocker mBlocker1;
    private final ReadStatusListener mListener;

    /**
     * This version operates both on mated and unmated results as determined
     * by the <code>blocker2</code> parameter. To be more specific it operates on the AS field
     * when supplied with a read blocker for each side, or on the XA field when only provided
     * with a single read blocker.
     *
     * @param blocker1 information about the min score for each read. (left only for unmated case).
     * @param listener each record written will set a flag in this listener
     */
    DummySamResultsFilter(MapQScoringReadBlocker blocker1, ReadStatusListener listener) {
      super(new MockSequencesReader(SequenceType.DNA, 21), new MockSequencesReader(SequenceType.DNA, 21), null, false, 0, false, false);
      mBlocker1 = blocker1;
      mListener = listener;
    }

    @Override
    protected String getName() {
      return "Dummy";
    }

    @Override
    protected SAMRecord filterRecord(SAMFileWriter samWriter, BinaryTempFileRecord rec, NamesInterface templateNames) {
      final int readId = rec.getReadId();
      final int flag = rec.getSamFlags() & 0xff;
      //assert (flag & SamBamRecord.SAM_READ_IS_PAIRED) == 0;
      assert (flag & SamBamConstants.SAM_READ_IS_UNMAPPED) == 0;
      final int score = rec.getAlignmentScore();
      if (score < 0) {
        // oops! no AS:score in this record
        throw new RuntimeException("SAM record has no " + SamUtils.ATTRIBUTE_ALIGNMENT_SCORE + ":score field!  readId=" + readId + " flag=" + flag);
      }
      //System.out.println("readId=" + readId + ", score=" + score);
      if (!mBlocker1.isBlocked1(readId, score)) {
        final int count = mBlocker1.getCount1(readId);
        final SAMRecord record = new SAMRecord(samWriter.getFileHeader());
        record.setReadName(String.valueOf(readId));
        record.setFlags(flag);
        record.setReferenceIndex(rec.getReferenceId());
        record.setAlignmentStart(rec.getStartPosition());
        record.setMappingQuality(255);
        record.setCigarString(new String(rec.getCigarString()));
        if (rec.isReadPaired()) {
          record.setMateReferenceIndex(record.getReferenceIndex());
          record.setMateAlignmentStart(rec.getMatePosition());
          record.setInferredInsertSize(rec.getTemplateLength());
        }
        addBaseAttributes(record, rec);
        record.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "G1");
        record.setAttribute(SamUtils.ATTRIBUTE_IH, count == 0 ? 1 : count);
        if (readId == 1) {
          record.setReadString("AGCT");
        } else if (readId == 9) {
          if (rec.getStartPosition() == 3) {
          record.setReadString("CTAG");
          } else {
            record.setReadString("TAGC");
          }
        } else {
          if (rec.getStartPosition() == 6) {
            record.setReadString("GCTA");
          } else {
            record.setReadString("AGCT");
          }
        }
//        + "1\t0\t0\t1\t20\t4M\t=\t0\t0\tAGCT\t<<<<\tNM:i:1\tAS:i:3\tRG:Z:G1\n"
//        + "9\t16\t0\t3\t30\t4M\t=\t0\t0\tCTAG\t<<<<\tMF:i:18\tAS:i:0\tRG:Z:G1\n"
//        + "9\t16\t0\t4\t30\t4M\t=\t0\t0\tTAGC\t<<<<\tMF:i:18\tAS:i:1\tRG:Z:G1\n"
//        + "21\t0\t0\t6\t30\t4M\t=\t0\t0\tGCTA\t<<<<\tMF:i:18\tAS:i:4\tRG:Z:G1\n"
//        + "21\t2\t0\t9\t30\t4M\t=\t0\t0\tAGCT\t<<<<\tMF:i:18\tAS:i:4\tRG:Z:G1\n";

        record.setBaseQualityString("<<<<");
        samWriter.addAlignment(record);
        mListener.addStatus(readId, ReadStatusTracker.UNMATED_FIRST);
        return record;
      }
      return null;
    }
  }
  protected static class StatusListener implements ReadStatusListener {
    private final int[] mStatus;

    public StatusListener(final int numRecords) {
      mStatus = new int[numRecords];
    }

    @Override
    public void addStatus(final int readId, final int status) {
      mStatus[readId] |= status;
    }

    public int getStatus(final int readId) {
      return mStatus[readId] & ReadStatusTracker.MAPPING_STATUS_MASK;
    }
  }

  public void testTimeCompute() {
    final DummySamResultsFilter filter = new DummySamResultsFilter(null, null);
    assertEquals(29L, filter.diffDivided(29726124));
    assertEquals(100L, filter.diffDivided(100000000L));
    assertEquals("4709.66", Utils.realFormat(filter.bytesPerSecond(29726124, 140), 2));
  }

  private static final String SAM_UNMATED_EXPECTED = ""
    + "1\t1\tchr20\t1\t255\t4=\t=\t0\t0\tAGCT\t<<<<\tAS:i:3\tNM:i:1\tRG:Z:G1\tIH:i:1\n"
    + "9\t17\tchr20\t3\t255\t4M\t=\t0\t0\tCTAG\t<<<<\tAS:i:0\tNM:i:0\tRG:Z:G1\tIH:i:2\n";


  private MockArraySequencesReader mTemplateReader;

  @Override
  protected void setUp() {
    mTemplateReader = new MockArraySequencesReader(SequenceType.DNA, new int[] {62435964}, new String[] {"chr20"})  ;
  }

  @Override
  protected void tearDown() {
    mTemplateReader = null;
  }

  private SAMFileHeader makeHeader() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMSequenceDictionary dict = header.getSequenceDictionary();
    dict.addSequence(new SAMSequenceRecord("chr20", 62435964));
    return header;
  }

  public void writeTempFile(File out) throws IOException {
    try (TempRecordWriter trw = new TempRecordWriterNio(FileUtils.createOutputStream(out, true))) {
      final BinaryTempFileRecord bar = new BinaryTempFileRecord(true, false, false, false);

      bar.setAlignmentScore(3);
      bar.setCigarString("4=".getBytes());
      bar.setMatePosition(0);
      bar.setMdString(new byte[0]);
      bar.setNumberMismatches(1);
      bar.setReadId(1);
      bar.setReferenceId(0);
      bar.setSamFlags((byte) 1);
      bar.setStartPosition(1);
      bar.setTemplateLength(0);
      trw.writeRecord(bar);

      bar.setAlignmentScore(0);
      bar.setCigarString("4M".getBytes());
      bar.setNumberMismatches(0);
      bar.setReadId(9);
      bar.setSamFlags((byte) 17);
      bar.setStartPosition(3);
      trw.writeRecord(bar);

      bar.setAlignmentScore(1);
      bar.setNumberMismatches(0);
      bar.setStartPosition(4);
      trw.writeRecord(bar);

      bar.setReadId(11);
      bar.setAlignmentScore(4);
      bar.setSamFlags((byte) 0);
      bar.setStartPosition(6);
      trw.writeRecord(bar);

      bar.setSamFlags((byte) 2);
      bar.setStartPosition(9);
      trw.writeRecord(bar);

      bar.setSentinelRecord();
      trw.writeRecord(bar);
    }
  }

  public void testFilterUnmated() throws IOException {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    try (PrintStream prLog = new PrintStream(log)) {
      Diagnostic.setLogStream(prLog);
      final int numReads = 100;
      final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 2);
      blocker.increment(1, 3);
      blocker.increment(11, 4);
      blocker.increment(11, 4);
      blocker.increment(11, 4); // read 21 is blocked for score=4
      blocker.increment(9, 1);
      blocker.increment(9, 1);
      blocker.increment(9, 1); // read 9 is blocked for score=1
      blocker.increment(9, 0);
      blocker.increment(9, 0); // read 9 is just not blocked for score=0
      try (final TestDirectory dir = new TestDirectory("unmatedSamFilter")) {
        final File in1 = File.createTempFile("sam", "_1.gz", dir);
        final File outFile = File.createTempFile("out", ".gz", dir);
        final StatusListener listener = new StatusListener(numReads);
        try (OutputStream out = new GZIPOutputStream(new FileOutputStream(outFile))) {

          writeTempFile(in1);
          final DummySamResultsFilter filter = new DummySamResultsFilter(blocker, listener);
          assertEquals("Dummy", filter.getName());
          filter.filterConcat(makeHeader(), out, null, null, mTemplateReader, false, in1);
        }
        final String contents = FileHelper.gzFileToString(outFile);
//        System.out.println("contents=" + contents);
        assertTrue(TestUtils.sameLines(SAM_UNMATED_EXPECTED, TestUtils.stripSAMHeader(contents), false));

        // now check that the listener has been updated correctly.
        for (int read = 0; read < numReads; ++read) {
          final int expect;
          switch (read) {
            case 1:
            case 9:
              expect = ReadStatusTracker.UNMATED_FIRST;
              break;
            default:
              expect = 0;
              break;
          }
          assertEquals("readId=" + read, expect, listener.getStatus(read));
        }
      }
    } finally {
      Diagnostic.setLogStream();
    }
    final String logString = log.toString();
    //System.err.println(logString);
    //String[] parts = logString.split("=");
    //    double len = Double.parseDouble(parts[7].substring(0, parts[7].indexOf("time")));
    //    double diff = Double.parseDouble(parts[8].substring(0, parts[8].indexOf("ms")));
    //    double result = Double.parseDouble(parts[9].substring(0, parts[8].indexOf(" ")));
    ////    System.err.println("result " + result);
    ////    System.err.println("diff " + ((len * 1.0e9 / (diff * 1000000)) + " diff + 1 " + (len * 1.0e9 / ((diff + 1) * 1000000))));
    //    double upper = len * 1.0e9 / (diff * 1000000);
    //    double lower = len * 1.0e9 / ((diff + 1) * 1000000);
    //    assertTrue(result < upper && result > lower);
    TestUtils.containsAll(logString, "Dummy SAM filter outputs 2/5 records",
      "filter concat file=", "bytes=", "time=", "ms", "bytes/sec=",
      "concat", "Dummy SAM filter outputs");
  }

  public void testFilterNoHeader() throws IOException {
    Diagnostic.setLogStream();
    final int numReads = 100;
    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 2);
    blocker.increment(1, 3);
    blocker.increment(11, 4);
    blocker.increment(11, 4);
    blocker.increment(11, 4); // read 11 is blocked for score=4
    blocker.increment(9, 1);
    blocker.increment(9, 1);
    blocker.increment(9, 1); // read 9 is blocked for score=1
    blocker.increment(9, 0);
    blocker.increment(9, 0); // read 9 is just not blocked for score=0
    try (final TestDirectory dir = new TestDirectory("unmatedSamFilter")) {
      final File in1 = File.createTempFile("sam", "_1.gz", dir);
      final File outFile = File.createTempFile("out", ".gz", dir);
      writeTempFile(in1);
      final StatusListener listener = new StatusListener(numReads);
      try (OutputStream out = new GZIPOutputStream(new FileOutputStream(outFile))) {

        final DummySamResultsFilter filter = new DummySamResultsFilter(blocker, listener);
        filter.setWriteHeader(false);
        assertEquals("Dummy", filter.getName());
        filter.filterConcat(makeHeader(), out, null, null, mTemplateReader, false, in1);
      }
      final String contents = FileHelper.gzFileToString(outFile);
      //System.out.println("contents=" + contents);
      assertTrue(TestUtils.sameLines(SAM_UNMATED_EXPECTED, contents, false));

      // now check that the listener has been updated correctly.
      for (int read = 0; read < numReads; ++read) {
        final int expect;
        switch (read) {
          case 1:
          case 9:
            expect = ReadStatusTracker.UNMATED_FIRST;
            break;
          default:
            expect = 0;
            break;
        }
        assertEquals("readId=" + read, expect, listener.getStatus(read));
      }
    }
  }

  public void testCalibrationCreation() throws IOException {
    CommandLine.clearCommandArgs();
    Diagnostic.setLogStream();
    final int numReads = 100;
    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 2);
    try (final TestDirectory dir = new TestDirectory()) {
      final File in1 = File.createTempFile("sam", "_1.gz", dir);
      final String templateString = ">chr20\nAGCTAGCTAGCTAGCTAGCT\n";
      //                                     12345678901234567890
      final File temp = new File(dir, "input");
      final File outFile = File.createTempFile("out", ".gz", dir);
      writeTempFile(in1);
      final File calFile = new File(outFile.getPath() + ".calibration");
      try (SequencesReader template = ReaderTestUtils.getReaderDNA(templateString, temp, new SdfId(0));
           OutputStream out = new GZIPOutputStream(new FileOutputStream(outFile));
           OutputStream calOut = new FileOutputStream(calFile)) {

        final StatusListener listener = new StatusListener(numReads);
        final DummySamResultsFilter filter = new DummySamResultsFilter(blocker, listener);
        assertEquals("Dummy", filter.getName());
        filter.filterConcat(makeHeader(), out, calOut, null, template, false, in1);
      }
      final String expected = "#CL\tnull" + StringUtils.LS
        + "@nh:G1\t0\t5" + StringUtils.LS
        + "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel" + StringUtils.LS
        + "G1\t27\tchr20\t20\t0\t0\t0" + StringUtils.LS;
      assertEquals(expected, CalibratorTest.stripVersion(FileUtils.fileToString(calFile)));
    }
  }

  public void testCalibrationCreationWithBedRegions() throws IOException {
    CommandLine.clearCommandArgs();
    Diagnostic.setLogStream();
    final int numReads = 100;
    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 2);
    try (final TestDirectory dir = new TestDirectory()) {
      final ReferenceRegions referenceRegions = BedUtils.regions(FileUtils.stringToFile("chr20\t2\t10", new File(dir, "regions.bed")));
      final File in1 = File.createTempFile("sam", "_1.gz", dir);
      final String templateString = ">chr20\nAGCTAGCTAGCTAGCTAGCT\n>chr21\nAAAAAACCCCTTTTTGGGGG\n";
      //                                     12345678901234567890          12345678901234567890
      final File temp = new File(dir, "input");
      final File outFile = File.createTempFile("out", ".gz", dir);
      writeTempFile(in1);
      final File calFile = new File(outFile.getPath() + ".calibration");
      try (SequencesReader template = ReaderTestUtils.getReaderDNA(templateString, temp, new SdfId(0));
        OutputStream out = new GZIPOutputStream(new FileOutputStream(outFile));
        OutputStream calOut = new FileOutputStream(calFile)) {

        final StatusListener listener = new StatusListener(numReads);
        final DummySamResultsFilter filter = new DummySamResultsFilter(blocker, listener);
        assertEquals("Dummy", filter.getName());
        filter.filterConcat(makeHeader(), out, calOut, referenceRegions, template, false, in1);
      }
      final String expected = "#CL\tnull" + StringUtils.LS
        + "@nh:G1\t0\t5" + StringUtils.LS
        + "@sequence\t8\tchr20" + StringUtils.LS
        + "@sequence\t0\tchr21" + StringUtils.LS
        + "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel" + StringUtils.LS
        + "G1\t27\tchr20\t16\t0\t0\t0" + StringUtils.LS;
      assertEquals(expected, CalibratorTest.stripVersion(FileUtils.fileToString(calFile)));
    }
  }
}
