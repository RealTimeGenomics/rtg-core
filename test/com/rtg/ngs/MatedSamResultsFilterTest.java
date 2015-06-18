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

import com.rtg.mode.SequenceType;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.TempRecordWriter;
import com.rtg.ngs.tempstage.TempRecordWriterNio;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import junit.framework.TestCase;

/**
 */
public class MatedSamResultsFilterTest extends TestCase {

  private static void writeRecords1(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, false, false);
    bpar.setReadId(1);
    bpar.setSamFlags((byte) 99);
    bpar.setReferenceId(0);
    bpar.setStartPosition(28833);
    bpar.setCigarString("4M5M".getBytes());
    bpar.setMatePosition(28993);
    bpar.setTemplateLength(195);
    bpar.setNumberMismatches(1);
    bpar.setComboScore(3);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
    bpar.setReadId(20);
    bpar.setSamFlags((byte) 147);
    bpar.setReferenceId(0);
    bpar.setStartPosition(28834);
    bpar.setCigarString("35M".getBytes());
    bpar.setMatePosition(28701);
    bpar.setTemplateLength(-168);
    bpar.setNumberMismatches(0);
    bpar.setComboScore(4);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }

  private static void writeRecords2(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, false, false);
    bpar.setReadId(30);
    bpar.setSamFlags((byte) 99);
    bpar.setReferenceId(0);
    bpar.setStartPosition(28833);
    bpar.setCigarString("4M5M".getBytes());
    bpar.setMatePosition(28993);
    bpar.setTemplateLength(195);
    bpar.setNumberMismatches(1);
    bpar.setComboScore(5);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
    bpar.setReadId(59);
    bpar.setSamFlags((byte) 147);
    bpar.setReferenceId(0);
    bpar.setStartPosition(28834);
    bpar.setCigarString("35M".getBytes());
    bpar.setMatePosition(28701);
    bpar.setTemplateLength(-168);
    bpar.setComboScore(6);
    bpar.setNumberMismatches(0);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }

  private static void writeRecords2Extra(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, false, false);
    bpar.setReadId(66);
    bpar.setSamFlags((byte) 147);
    bpar.setReferenceId(0);
    bpar.setStartPosition(28834);
    bpar.setCigarString("35M".getBytes());
    bpar.setMatePosition(28701);
    bpar.setTemplateLength(-168);
    bpar.setComboScore(6);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
    bpar.setReadId(67);
    bpar.setSamFlags((byte) 147);
    bpar.setReferenceId(0);
    bpar.setStartPosition(28834);
    bpar.setCigarString("35M".getBytes());
    bpar.setMatePosition(28701);
    bpar.setTemplateLength(-168);
    bpar.setComboScore(6);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }

  private static void writeSentinel(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, false, false);
    bpar.setSentinelRecord();
    trw.writeRecord(bpar);
  }

  private static final String EXPECTED = ""
    + "1\t99\tchr20\t28833\t55\t4M5M\t=\t28993\t195\tACGT\t<<<<\tAS:i:0\tNM:i:1\tXA:i:3\tIH:i:1\tNH:i:1\n"
    + "20\t403\tchr20\t28834\t0\t35M\t=\t28701\t-168\tACGT\t<<<<\tAS:i:0\tNM:i:0\tXA:i:4\tIH:i:4\tNH:i:4\n"
    + "59\t147\tchr20\t28834\t55\t35M\t=\t28701\t-168\tACGT\t<<<<\tAS:i:0\tNM:i:0\tXA:i:6\tIH:i:1\tNH:i:1\n";

  private MockArraySequencesReader mTemplateReader;

  private static MockSequencesReader getReader() {
    return new MockSequencesReader(SequenceType.DNA, 100) {
      @Override
      public int read(long index, byte[] out) {
        return read(index, out, 0, 4);
      }

      @Override
      public int read(long index, byte[] out, int start, int length) {
        System.arraycopy(new byte[] {(byte) 1, (byte) 2, (byte) 3, (byte) 4}, start, out, 0, length);
        return 4;
      }

      @Override
      public boolean hasQualityData() {
        return true;
      }

      @Override
      public int readQuality(long sequenceIndex, byte[] dest) {
        return readQuality(sequenceIndex, dest, 0, 4);
      }

      @Override
      public int readQuality(long sequenceIndex, byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException {
        System.arraycopy(new byte[] {(byte) ('<' - '!'), (byte) ('<' - '!'), (byte) ('<' - '!'), (byte) ('<' - '!')}, start, dest, 0, length);
        return 4;
      }
    };
  }


  @Override
  protected void setUp() {
    mTemplateReader = new MockArraySequencesReader(SequenceType.DNA, new int[] {62435964}, new String[] {"chr20"})  ;
  }

  @Override
  protected void tearDown() {
    mTemplateReader = null;
  }

  public void testFilterConcatGzipped() throws IOException {
    runFilterConcat(true, false);
  }

  public void testFilterConcat() throws IOException {
    runFilterConcat(false, false);
  }

  public void testFilterConcatGzippedDeleteIntermediate() throws IOException {
    runFilterConcat(true, true);
  }

  public void testFilterConcatDeleteIntermediate() throws IOException {
    runFilterConcat(false, true);
  }

  public void runFilterConcat(final boolean gzipped, boolean deleteIntermediate) throws IOException {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    try (PrintStream printLog = new PrintStream(log)) {
      Diagnostic.setLogStream(printLog);
      final int numReads = 100;
      final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 4);
      blocker.increment(1, 3);
      blocker.increment(20, 4);
      blocker.increment(20, 4);
      blocker.increment(20, 4);
      blocker.increment(20, 4);
      blocker.increment(30, 5);
      blocker.increment(30, 5);
      blocker.increment(30, 5);
      blocker.increment(30, 5);
      blocker.increment(30, 5);
      blocker.increment(59, 6);
      final File dir = FileUtils.createTempDir("test", "matedSamFilter");
      OutputStream out = null;
      try {
        final File in1 = File.createTempFile("sam", "_1.gz", dir);

        try (TempRecordWriter trw = new TempRecordWriterNio(FileUtils.createOutputStream(in1, true))) {
          writeRecords1(trw);
          writeSentinel(trw);
        }
        final File in2 = File.createTempFile("sam", "_2.gz", dir);
        try (TempRecordWriter trw2 = new TempRecordWriterNio(FileUtils.createOutputStream(in2, true))) {
          writeRecords2(trw2);
          writeRecords2Extra(trw2);
          writeSentinel(trw2);
        }
        //FileHelper.stringToGzFile(SAM1, in1);
        //FileHelper.stringToGzFile(SAM2 + SAM2_EXTRA, in2);
        final File outFile = File.createTempFile("out", ".gz", dir);
        out = new FileOutputStream(outFile);
        if (gzipped) {
          out = new GZIPOutputStream(out);
        }
        final ReadBlocker freqBlockerLeft = new ReadBlocker(numReads, 2);
        final ReadBlocker freqBlockerRight = new ReadBlocker(numReads, 2);
        freqBlockerLeft.increment(66);
        freqBlockerLeft.increment(66);
        freqBlockerRight.increment(67);
        freqBlockerRight.increment(67);
        final MatedSamResultsFilter filter = new MatedSamResultsFilter(blocker, freqBlockerLeft, freqBlockerRight, getReader(), getReader(), false, 0, null, false);
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addSequence(new SAMSequenceRecord("chr20", 62435964));
        filter.setHeader(header);
        filter.filterConcat(header, out, null, null, mTemplateReader, deleteIntermediate, in1, in2);

        if (!deleteIntermediate) {
          assertTrue(in1.exists());
          assertTrue(in2.exists());
        } else {
          assertFalse(in1.exists());
          assertFalse(in2.exists());
        }
        out.close();
        final String contents = gzipped
          ? FileHelper.gzFileToString(outFile)
          : FileUtils.fileToString(outFile);
        //System.out.println("contents=" + contents);
        assertTrue(TestUtils.sameLines(EXPECTED, TestUtils.stripSAMHeader(contents), false));

      } finally {
        if (out != null) {
          out.close();
        }
        assertTrue(FileHelper.deleteAll(dir));
        Diagnostic.closeLog();
        Diagnostic.setLogStream();
        //System.out.println(log.toString());
      }
    }
    assertTrue(log.toString().contains("Mated SAM filter outputs 3/6 records"));
  }

  public void testFilterConcatEmpty() throws IOException {
    final MemoryPrintStream input = new MemoryPrintStream();
    try (TempRecordWriter trw = new TempRecordWriterNio(new GZIPOutputStream(input.outputStream()))) {
      writeSentinel(trw);
    }
    final OutputStream log = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(log));
    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(100, 4);
    final File dir = FileUtils.createTempDir("test", "matedSamFilter");
    OutputStream out = null;
    try {
      final File in1 = File.createTempFile("sam", "_1.gz", dir);
      final File in2 = File.createTempFile("sam", "_2.gz", dir);
      FileUtils.byteArrayToFile(input.toByteArray(), in1);
      FileUtils.byteArrayToFile(input.toByteArray(), in2);
      final File outFile = File.createTempFile("out", ".gz", dir);
      out = new GZIPOutputStream(new FileOutputStream(outFile));

      final MatedSamResultsFilter filter = new MatedSamResultsFilter(blocker, new ReadBlocker(100, 255), new ReadBlocker(100, 255), null, null, false, 0, null, false);
      final SAMFileHeader header = new SAMFileHeader();
      header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      header.addSequence(new SAMSequenceRecord("chr20", 62435964));
      filter.filterConcat(header, out, null, null, mTemplateReader, false, in1, in2);
      out.close();
      final String contents = FileHelper.gzFileToString(outFile);
      //System.out.println("contents=" + contents);
      assertTrue(contents.startsWith("@HD"));
      assertTrue(TestUtils.sameLines("", TestUtils.stripSAMHeader(contents), false));
    } finally {
      if (out != null) {
        out.close();
      }
      Diagnostic.setLogStream();
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  /**
   * A main method that can be used to measure the speed of the filtering
   * on large SAM files (up to 1 million reads)
   *
   * @param args infile.gz outfile.gz
   * @throws IOException on IO error
   */
  public static void main(final String[] args) throws IOException {
    if (args.length != 2) {
      System.out.println("Arguments in.sam.gz out.sam.gz");
      System.exit(1);
    }
    final int numReads = 1000000;
    final File in1 = new File(args[0]);
    final File out = new File(args[1]);
    try (OutputStream outStream = new GZIPOutputStream(new FileOutputStream(out))) {
      final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 4);
      final MatedSamResultsFilter filter = new MatedSamResultsFilter(blocker, new ReadBlocker(numReads, 255), new ReadBlocker(numReads, 255), null, null, false, 0, null, false);
      filter.filterConcat(null, outStream, null, null, null, false, in1);
    }
  }

  private void checkCG(final byte[] recordsIn, final String samOut, String readGroupId) throws Exception {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    try (PrintStream pr = new PrintStream(ba)) {
      Diagnostic.setLogStream(pr);
      try {
        final File mainOut = FileUtils.createTempDir("cgmap", "test");
        try {
          final File cgleft = new File(mainOut, "cgleft");
          final File cgright = new File(mainOut, "cgright");
          try (SequencesReader lr = ReaderTestUtils.getReaderDNAFastqCG("@s1\nctggttaaaatatgaagtgaccaccatgcttgaga\n+\nABCDEFGHIJKLMNOPQRSTUVWXYZ123456789\n", cgleft, PrereadArm.LEFT); SequencesReader rr = ReaderTestUtils.getReaderDNAFastqCG("@s1\ntctcaagcatggtggtcacttcatattttaaccag\n+\nABCDEFGHIJKLMNOPQRSTUVWXYZ123456789\n", cgright, PrereadArm.RIGHT)) {
            final int numReads = 100;
            final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 4);
            final File in1 = File.createTempFile("sam", "_1.gz", mainOut);
            FileUtils.byteArrayToFile(recordsIn, in1);
            final File outFile = File.createTempFile("out", ".gz", mainOut);
            try (OutputStream out = new FileOutputStream(outFile)) {
              for (int k = 0; k < numReads; k++) {
                blocker.increment(k, 50);
              }
              final MatedSamResultsFilter filter = new MatedSamResultsFilter(blocker, new ReadBlocker(numReads, 255), new ReadBlocker(numReads, 255), lr, rr, true, 0, readGroupId, false);
              final SAMFileHeader header = new SAMFileHeader();
              header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
              header.addSequence(new SAMSequenceRecord("chr20", 62435964));
              filter.setHeader(header); //XXX this is dumb V
              filter.filterConcat(header, out, null, null, mTemplateReader, false, in1);
            }
            final String contents = FileUtils.fileToString(outFile);
            // System.out.println("contents=" + contents);
            assertTrue(TestUtils.sameLines(samOut, TestUtils.stripSAMHeader(contents), false));
          }
        } finally {
          FileHelper.deleteAll(mainOut);
        }
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

  private static void writeCgRecords1(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, true, false);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 67);
    bpar.setReferenceId(0);
    bpar.setStartPosition(1);
    bpar.setCigarString("24=5N2X8=".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(2);
    bpar.setCgReadString("CTGGTAAAATATGAAGTGACCACCATGCTTGAGA".getBytes());
    bpar.setReadDeltaString("AT".getBytes());
    bpar.setSuperCigarString("5=1B20=5N2X8=".getBytes());
    bpar.setComboScore(3);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 147);
    bpar.setReferenceId(0);
    bpar.setStartPosition(1);
    bpar.setCigarString("24=5N2X8=".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(2);
    bpar.setCgReadString("CTGGTAAAATATGAAGTGACCACCATGCTTGAGA".getBytes());
    bpar.setReadDeltaString("AT".getBytes());
    bpar.setSuperCigarString("5=1B20=5N2X8=".getBytes());
    bpar.setComboScore(3);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }
  private static final String EXPECTED_CG_SUPERCIGAR_1 = ""
    + "0\t67\tchr20\t1\t55\t24=5N2X8=\t=\t0\t0\tCTGGTAAAATATGAAGTGACCACCATGCTTGAGA\tABCDEGHIJKLMNOPQRSTUVWXYZ123456789\tAS:i:2\tXU:Z:5=1B20=5N2X8=\tXR:Z:AT\tXA:i:3\tXQ:Z:F\tIH:i:1\tNH:i:1\n"
    + "0\t147\tchr20\t1\t55\t24=5N2X8=\t=\t0\t0\tCTGGTAAAATATGAAGTGACCACCATGCTTGAGA\t98765321ZYXWVUTSRQPONMLKJIHGFEDCBA\tAS:i:2\tXU:Z:5=1B20=5N2X8=\tXR:Z:AT\tXA:i:3\tXQ:Z:4\tIH:i:1\tNH:i:1\n";

  /** test that quality and XQ (extra quality) fields are calculated correctly from super cigar. */
  public void testCGSuperCigar() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    try (TempRecordWriter trw = new TempRecordWriterNio(new GZIPOutputStream(mps.outputStream()))) {
      writeCgRecords1(trw);
      writeSentinel(trw);
    }
    checkCG(mps.toByteArray(), EXPECTED_CG_SUPERCIGAR_1, null);
  }

  private static void writeCgRecords2(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, true, false);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 67);
    bpar.setReferenceId(0);
    bpar.setStartPosition(1);
    bpar.setCigarString("blah".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(2);
    bpar.setCgReadString("CTGGTAAAATATGAAGTGACCACCATGCTTGA".getBytes());
    bpar.setReadDeltaString("ATT".getBytes());
    bpar.setSuperCigarString("2=1X2=3B3=2I25=".getBytes());
    bpar.setComboScore(3);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 131);
    bpar.setReferenceId(0);
    bpar.setStartPosition(1);
    bpar.setCigarString("boo".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(2);
    bpar.setCgReadString("CTGGTAAAATATGAAGTGACCACCATGCTTG".getBytes());
    bpar.setReadDeltaString("ACCGGG".getBytes());
    bpar.setSuperCigarString("2=1X4=2I14=3I4=4B2=1D3=".getBytes());
    bpar.setComboScore(3);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }


  private static final String EXPECTED_CG_SUPERCIGAR_2 = ""
    + "0\t67\tchr20\t1\t55\tblah\t=\t0\t0\tCTGGTAAAATATGAAGTGACCACCATGCTTGA\tABCDEIJKLMNOPQRSTUVWXYZ123456789\tAS:i:2\tXU:Z:2=1X2=3B3=2I25=\tXR:Z:ATT\tXA:i:3\tXQ:Z:FGH\tIH:i:1\tNH:i:1\n"
    + "0\t131\tchr20\t1\t55\tboo\t=\t0\t0\tCTGGTAAAATATGAAGTGACCACCATGCTTG\tABCDEFGHIJKLMNOPQRSTUVWXYZ56789\tAS:i:2\tXU:Z:2=1X4=2I14=3I4=4B2=1D3=\tXR:Z:ACCGGG\tXA:i:3\tXQ:Z:1234\tIH:i:1\tNH:i:1\n";

  /** test that quality and XQ (extra quality) fields are calculated correctly from more complex super cigars. */
  public void testCGSuperCigar2() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    try (TempRecordWriter trw = new TempRecordWriterNio(new GZIPOutputStream(mps.outputStream()))) {
      writeCgRecords2(trw);
      writeSentinel(trw);
    }
    checkCG(mps.toByteArray(), EXPECTED_CG_SUPERCIGAR_2, null);
  }

  private static void writeCgRecords3(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, true, false);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 67);
    bpar.setReferenceId(0);
    bpar.setStartPosition(1);
    bpar.setCigarString("blah".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(2);
    bpar.setCgReadString("CTGGTANAAATATGAAGTGACCACCATGCTTGAGA".getBytes());
    bpar.setReadDeltaString("TGGN".getBytes());
    bpar.setSuperCigarString("1=3I1=3D3B1=1I28=".getBytes());
    bpar.setComboScore(4);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 131);
    bpar.setReferenceId(0);
    bpar.setStartPosition(1);
    bpar.setCigarString("boo".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(2);
    bpar.setCgReadString("CTGGTAAAAATATGAAGTGACCACCATGCTTGAGA".getBytes());
    bpar.setReadDeltaString(new byte[0]);
    bpar.setSuperCigarString(new byte[0]);
    bpar.setComboScore(4);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }

  private static final String EXPECTED_CG_SUPERCIGAR_3 = ""
    + "0\t67\tchr20\t1\t55\tblah\t=\t0\t0\tCTGGTANAAATATGAAGTGACCACCATGCTTGAGA\tABCDEFGHIJKLMNOPQRSTUVWXYZ123456789\tAS:i:2\tXU:Z:1=3I1=3D3B1=1I28=\tXR:Z:TGGN\tXA:i:4\tIH:i:1\tNH:i:1\n"
    + "0\t131\tchr20\t1\t55\tboo\t=\t0\t0\tCTGGTAAAAATATGAAGTGACCACCATGCTTGAGA\tABCDEFGHIJKLMNOPQRSTUVWXYZ123456789\tAS:i:2\tXA:i:4\tIH:i:1\tNH:i:1\n";

  /** test that quality field is added, but no XQ is added when no CG overlap (so no super cigar). */
  public void testCGSuperCigar3() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    try (TempRecordWriter trw = new TempRecordWriterNio(new GZIPOutputStream(mps.outputStream()))) {
      writeCgRecords3(trw);
      writeSentinel(trw);
    }
    checkCG(mps.toByteArray(), EXPECTED_CG_SUPERCIGAR_3, null);
  }

  public void testInvalidQualityLength() throws Exception {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, true, false);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 131);
    bpar.setReferenceId(0);
    bpar.setStartPosition(17617);
    bpar.setCigarString("10=6N18=2I1=".getBytes());
    bpar.setMatePosition(0);
    bpar.setTemplateLength(0);
    bpar.setAlignmentScore(3);
    bpar.setCgReadString("CGAGTATGGACCTCATAAGTCGGCCCCTGGT".getBytes());
    bpar.setReadDeltaString("GG".getBytes());
    bpar.setSuperCigarString("10=6N20=4B2=2I1=".getBytes());
    bpar.setComboScore(4);
    bpar.setMdString(new byte[0]);
    final MemoryPrintStream input = new MemoryPrintStream();
    try (TempRecordWriter trw = new TempRecordWriterNio(new GZIPOutputStream(input.outputStream()))) {
      trw.writeRecord(bpar);
      writeSentinel(trw);
    }
    final String exp = ""
      + "0\t131\tchr20\t17617\t55\t10=6N18=2I1=\t=\t0\t0\tCGAGTATGGACCTCATAAGTCGGCCCCTGGT\tUZURRSWZSZRYSQRWSQSSURSWSWTSUQV\tAS:i:3\tXU:Z:10=6N20=4B2=2I1=\tXR:Z:GG\tXA:i:4\tXQ:Z:UXTQ\tIH:i:1\tNH:i:1\n";

    final File mainOut = FileUtils.createTempDir("cgmap", "test");
    try {
      final File cgleft = new File(mainOut, "cgleft");
      final File cgright = new File(mainOut, "cgright");
      try (SequencesReader lr = ReaderTestUtils.getReaderDNAFastqCG("@s1\nctggttaaaatatgaagtgaccaccatgcttgaga\n+\nABCDEFGHIJKLMNOPQRSTUVWXYZ123456789\n", cgleft, PrereadArm.LEFT); SequencesReader rr = ReaderTestUtils.getReaderDNAFastqCG("@s1\ncgagtatggacctcataagtcggccccttcctggt\n+\nUZURRSWZSZRYSQRWSQSSURSWSWUXTQTSUQV\n", cgright, PrereadArm.RIGHT)) {
        final int numReads = 100;
        final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 4);
        final File in1 = File.createTempFile("sam", "_1.gz", mainOut);
        FileUtils.byteArrayToFile(input.toByteArray(), in1);
        final File outFile = File.createTempFile("out", ".gz", mainOut);
        try (OutputStream out = new FileOutputStream(outFile)) {
          for (int k = 0; k < numReads; k++) {
            blocker.increment(k, 50);
          }
          final MatedSamResultsFilter filter = new MatedSamResultsFilter(blocker, new ReadBlocker(numReads, 255), new ReadBlocker(numReads, 255), lr, rr, true, 0, null, false);
          final SAMFileHeader header = new SAMFileHeader();
          header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
          header.addSequence(new SAMSequenceRecord("chr20", 62435964));
          filter.setHeader(header);
          filter.filterConcat(header, out, null, null, mTemplateReader, false, in1);

          final String contents = FileUtils.fileToString(outFile);
          assertTrue(TestUtils.sameLines(exp, TestUtils.stripSAMHeader(contents), false));
        }
      }
    } finally {
      FileHelper.deleteAll(mainOut);
    }
  }

  private static void writeCgRecordsBug1295(TempRecordWriter trw) throws IOException {
    final BinaryTempFileRecord bpar = new BinaryTempFileRecord(true, false, true, false);
    bpar.setReadId(0);
    bpar.setSamFlags((byte) 179);
    bpar.setReferenceId(0);
    bpar.setStartPosition(782703);
    bpar.setCigarString("3=2I16=6N10=".getBytes());
    bpar.setMatePosition(783157);
    bpar.setTemplateLength(455);
    bpar.setAlignmentScore(3);
    bpar.setCgReadString("TGGCCAGAAAATGCAGAACAAAAGAACAGGG".getBytes());
    bpar.setReadDeltaString("CC".getBytes());
    bpar.setSuperCigarString("3=2I4B20=6N10=".getBytes());
    bpar.setComboScore(3);
    bpar.setMdString(new byte[0]);
    trw.writeRecord(bpar);
  }
  private static final String EXPECTED_CG_SUPERCIGAR_BUG1295 = ""
    + "0 179 chr20 782703  55  3=2I16=6N10=  = 783157  455 TGGCCAGAAAATGCAGAACAAAAGAACAGGG 98765ZYXWVUTSRQPONMLKJIHGFEDCBA  AS:i:3  XU:Z:3=2I4B20=6N10= XR:Z:CC XA:i:3  RG:Z:GS000005305  XQ:Z:4321  IH:i:1  NH:i:1\n".replaceAll(" +", "\t")
    ;

  /** test that quality and XQ (extra quality) fields are calculated correctly from more complex super cigars. */
  public void testCGQualLengthBug1295() throws Exception {
    final MemoryPrintStream input = new MemoryPrintStream();
    try (TempRecordWriter trw = new TempRecordWriterNio(new GZIPOutputStream(input.outputStream()))) {
      writeCgRecordsBug1295(trw);
      writeSentinel(trw);
    }
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      checkCG(input.toByteArray(), EXPECTED_CG_SUPERCIGAR_BUG1295, "GS000005305");
    } catch (final NoTalkbackSlimException ex) {
      System.err.println("ex: " + mps.toString());
      throw ex;
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
