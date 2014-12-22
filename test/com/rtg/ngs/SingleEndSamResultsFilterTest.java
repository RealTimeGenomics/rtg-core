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
import com.rtg.ngs.DummySamResultsFilterTest.StatusListener;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.TempRecordWriter;
import com.rtg.ngs.tempstage.TempRecordWriterNio;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.PrereadType;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class SingleEndSamResultsFilterTest extends TestCase {

  private static final String SAM_UNMATED_EXPECTED = "@HD\tVN:1.4\tSO:unsorted\n" + "@SQ\tSN:chr20\tLN:62435964\n"
  + "3\t272\tchr20\t28734\t0\t35M\t*\t0\t0\tACCT\t<<<<\tAS:i:0\tNM:i:0\tIH:i:2\tNH:i:2\n"
  + "1\t0\tchr20\t28833\t37\t4M5M\t*\t0\t0\tAGCT\t<<<<\tAS:i:3\tNM:i:1\tIH:i:1\tNH:i:1\n";


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
    final TempRecordWriter trw = new TempRecordWriterNio(FileUtils.createOutputStream(out, true, false));
    try {
      final BinaryTempFileRecord bar = new BinaryTempFileRecord(false, false, false, false);

//      + "3\t16\t0\t28734\t30\t35M\t=\t0\t0\tACCT\t<<<<\tAS:i:0\tNM:i:1\n"
      bar.setAlignmentScore(0);
      bar.setCigarString("35M".getBytes());
      bar.setMdString(new byte[0]);
      bar.setNumberMismatches(0);
      bar.setReadId(3);
      bar.setReferenceId(0);
      bar.setSamFlags((byte) 16);
      bar.setStartPosition(28734);
      trw.writeRecord(bar);

//      + "1\t0\t0\t28833\t20\t4M5M\t=\t0\t0\tAGCT\t<<<<\tNM:i:1\tAS:i:3\n"
      bar.setAlignmentScore(3);
      bar.setCigarString("4M5M".getBytes());
      bar.setNumberMismatches(1);
      bar.setReadId(1);
      bar.setSamFlags((byte) 0);
      bar.setStartPosition(28833);
      trw.writeRecord(bar);

//      + "3\t16\t0\t28834\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:1\n"
      bar.setAlignmentScore(1);
      bar.setCigarString("35M".getBytes());
      bar.setNumberMismatches(0);
      bar.setReadId(3);
      bar.setSamFlags((byte) 16);
      bar.setStartPosition(28834);
      trw.writeRecord(bar);

//      + "20\t0\t0\t28934\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:4\n"
      bar.setAlignmentScore(4);
      bar.setReadId(20);
      bar.setSamFlags((byte) 0);
      bar.setStartPosition(28934);
      trw.writeRecord(bar);

//      + "66\t0\t0\t28934\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:4\n";
      bar.setReadId(66);
      trw.writeRecord(bar);

      bar.setSentinelRecord();
      trw.writeRecord(bar);
    } finally {
      trw.close();
    }
  }

  public void testFilterUnmated() throws IOException {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    try (PrintStream prLog = new PrintStream(log)) {
      Diagnostic.setLogStream(prLog);
      final int numReads = 100;
      final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker(numReads, 2);
      blocker.increment(1, 3);
      blocker.increment(20, 4);
      blocker.increment(20, 4);
      blocker.increment(20, 4); // read 20 is blocked for score=4
      blocker.increment(3, 1);
      blocker.increment(3, 1);
      blocker.increment(3, 1); // read 3 is blocked for score=1
      blocker.increment(3, 0);
      blocker.increment(3, 0); // read 3 is just not blocked for score=0
      final File dir = FileUtils.createTempDir("test", "unmatedSamFilter");
      OutputStream out = null;
      try {
        final File in1 = File.createTempFile("sam", "_1.gz", dir);
        writeTempFile(in1);
        final File outFile = File.createTempFile("out", ".gz", dir);
        out = new GZIPOutputStream(new FileOutputStream(outFile));

        final StatusListener listener = new StatusListener(numReads);
        final ReadBlocker freqBlocker = new ReadBlocker(numReads, 2);
        freqBlocker.increment(66);
        freqBlocker.increment(66);

        final MockSequencesReader msr = new MockSequencesReader(SequenceType.DNA) {
          @Override
          public PrereadType getPrereadType() {
            return PrereadType.UNKNOWN;
          }

          @Override
          public boolean hasQualityData() {
            return true;
          }

          @Override
          public int read(long sequenceIndex, byte[] dataOut) {
            dataOut[0] = 1;
            if (sequenceIndex == 3) {
              dataOut[1] = 3;
              dataOut[2] = 3;
            } else {
              dataOut[1] = 3;
              dataOut[2] = 2;
            }
            dataOut[3] = 4;
            return 4;
          }

          @Override
          public int readQuality(final long sequenceIndex, final byte[] dest) {
            dest[0] = dest[1] = dest[2] = dest[3] = '<' - 33;
            return 4;
          }
        };

        final SingleEndSamResultsFilter filter = new SingleEndSamResultsFilter(blocker, freqBlocker, listener, 0, msr, null, false);
        assertEquals("Alignment", filter.getName());
        filter.filterConcat(makeHeader(), out, null, null, mTemplateReader, false, in1);
        out.close();
        final String contents = FileHelper.gzFileToString(outFile);
        //System.out.println("contents=" + contents);
        assertTrue(TestUtils.sameLines(SAM_UNMATED_EXPECTED, contents, false));

        // now check that the listener has been updated correctly.
        for (int read = 0; read < numReads; read++) {
          final int expect;
          switch (read) {
            case 1:
            case 3:
              expect = ReadStatusTracker.UNMATED_FIRST;
              break;
            default:
              expect = 0;
              break;
          }
          assertEquals("readId=" + read, expect, listener.getStatus(read));
        }
      } finally {
        if (out != null) {
          out.close();
        }
        assertTrue(FileHelper.deleteAll(dir));
      }
    }
    final String logString = log.toString();
    //System.err.println(logString);
    TestUtils.containsAll(logString, "Alignment SAM filter outputs 2/5 records");
  }
}
