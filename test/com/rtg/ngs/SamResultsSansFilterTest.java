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

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.SequenceType;
import com.rtg.ngs.DummySamResultsFilterTest.StatusListener;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.TempRecordWriter;
import com.rtg.ngs.tempstage.TempRecordWriterNio;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;


/**
 */
public class SamResultsSansFilterTest extends AbstractNanoTest {

  public void writeTempFile(File out) throws IOException {
    try (TempRecordWriter trw = new TempRecordWriterNio(FileUtils.createOutputStream(out, true, false))) {
      final BinaryTempFileRecord bar = new BinaryTempFileRecord(true, false, false, true);

//  + "1\t65\t0\t28833\t20\t4M5M\t=\t0\t0\tAGCT\t<<<<\tNM:i:1\tAS:i:3\n"
      bar.setAlignmentScore(3);
      bar.setCigarString("4M5M".getBytes());
      bar.setNumberMismatches(1);
      bar.setReadId(1);
      bar.setSamFlags((byte) 65);
      bar.setStartPosition(28833);
      bar.setReferenceId(0);
      bar.setMdString(new byte[0]);
      bar.setMatePosition(0);
      trw.writeRecord(bar);

//  + "3\t129\t0\t28734\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:0\n"
      bar.setAlignmentScore(0);
      bar.setCigarString("35M".getBytes());
      bar.setNumberMismatches(0);
      bar.setReadId(3);
      bar.setSamFlags((byte) 129);
      bar.setStartPosition(28734);
      trw.writeRecord(bar);

//  + "3\t129\t0\t28834\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:1\n"
      bar.setAlignmentScore(1);
      bar.setNumberMismatches(0);
      bar.setReadId(3);
      bar.setStartPosition(28834);
      trw.writeRecord(bar);

//  + "20\t65\t0\t28934\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:4\n"
      bar.setAlignmentScore(4);
      bar.setReadId(20);
      bar.setSamFlags((byte) 65);
      bar.setStartPosition(28934);
      trw.writeRecord(bar);

//  + "66\t65\t0\t28934\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:4\n"
      bar.setReadId(66);
      trw.writeRecord(bar);

//  + "67\t129\t0\t28934\t30\t35M\t=\t0\t0\tACCT\t<<<<\tMF:i:18\tAS:i:4\n";
      bar.setReadId(67);
      bar.setSamFlags((byte) 129);
      trw.writeRecord(bar);

      bar.setSentinelRecord();
      trw.writeRecord(bar);
    }
  }


  private MockArraySequencesReader mTemplateReader;

  @Override
  public void setUp() throws IOException {
    mTemplateReader = new MockArraySequencesReader(SequenceType.DNA, new int[] {62435964}, new String[] {"chr20"});
    super.setUp();
  }

  @Override
  public void tearDown() throws IOException {
    mTemplateReader = null;
    super.tearDown();
  }

  public void testFilterUnmated() throws IOException {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    try (PrintStream prLog = new PrintStream(log)) {
      Diagnostic.setLogStream(prLog);
      final int numReads = 100;
      final File dir = FileUtils.createTempDir("test", "unmatedSamFilter");
      OutputStream out = null;
      try {
        final File in1 = File.createTempFile("sam", "_1.gz", dir);
        writeTempFile(in1);
        final File outFile = File.createTempFile("out", ".gz", dir);
        out = new GZIPOutputStream(new FileOutputStream(outFile));

        final StatusListener listener = new StatusListener(numReads);

        try (SequencesReader mr = new MockSequencesReader(SequenceType.DNA, 80, 4) {
          @Override
          public int read(long index, byte[] out) {
            out[0] = 1;
            if (index == 1) {
              out[1] = 3;
            } else {
              out[1] = 2;
            }
            out[2] = 2;
            out[3] = 4;
            return 4;
          }

          @Override
          public PrereadType getPrereadType() {
            return PrereadType.UNKNOWN;
          }

          @Override
          public boolean hasQualityData() {
            return true;
          }

          @Override
          public int readQuality(long sequenceIndex, byte[] dest) {
            dest[0] = dest[1] = dest[2] = dest[3] = '<' - 33;
            return 4;
          }
        }) {
          final SamResultsSansFilter filter = new SamResultsSansFilter(listener, 0, mr, mr, null, false);
          assertEquals("Sans", filter.getName());

          filter.filterConcat(makeHeader(), out, null, null, mTemplateReader, false, in1);
          out.close();
          final String contents = FileHelper.gzFileToString(outFile);
          mNano.check("srsf-unmated", contents, false);

          // now check that the listener has been updated correctly.
          for (int read = 0; read < numReads; read++) {
            final int expect;
            switch (read) {
              case 1:
              case 20:
              case 66:
                expect = ReadStatusTracker.UNMATED_FIRST;
                break;
              case 3:
              case 67:
                expect = ReadStatusTracker.UNMATED_SECOND;
                break;
              default:
                expect = 0;
                break;
            }
            assertEquals("readId=" + read, expect, listener.getStatus(read));
          }
        }
      } finally {
        if (out != null) {
          out.close();
        }
        assertTrue(FileHelper.deleteAll(dir));
      }
      Diagnostic.setLogStream(prLog);
    }
    final String logString = log.toString();
    //System.err.println(logString);
    TestUtils.containsAll(logString, "Sans SAM filter outputs 6/6 records");
  }

  private SAMFileHeader makeHeader() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMSequenceDictionary dict = header.getSequenceDictionary();
    dict.addSequence(new SAMSequenceRecord("chr20", 62435964));
    return header;
  }
}
