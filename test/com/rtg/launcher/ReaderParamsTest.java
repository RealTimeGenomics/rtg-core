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
package com.rtg.launcher;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesWriter;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ReaderParamsTest extends TestCase {

  private File mDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final String SEQ_DNA_A = ""
    + ">x" + StringUtils.LS
    + "actg" + StringUtils.LS
    + ">y" + StringUtils.LS
    + "actg" + StringUtils.LS
    + ">z" + StringUtils.LS
    + "actg" + StringUtils.LS
    + ">a" + StringUtils.LS
    + "actg" + StringUtils.LS;

  ReaderParams getParams()  throws IOException {
    final File dir = getDir();
    return DefaultReaderParamsTest.createDefaultReaderParams(dir, LongRange.NONE, SequenceMode.BIDIRECTIONAL);
  }

  private File getDir() throws IOException {
    final File dir = FileHelper.createTempDirectory(mDir);
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(SEQ_DNA_A.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams,
        new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, false);
    sequenceWriter.processSequences();
    return dir;
  }

  public void testEquals() throws IOException, ClassNotFoundException {
    final File dira = getDir();
    final File dirb = getDir();
    final ReaderParams a1 = DefaultReaderParamsTest.createDefaultReaderParams(dira, LongRange.NONE, SequenceMode.BIDIRECTIONAL);
    final ReaderParams a2 = DefaultReaderParamsTest.createDefaultReaderParams(dira, LongRange.NONE, SequenceMode.BIDIRECTIONAL);
    final ReaderParams b1 = DefaultReaderParamsTest.createDefaultReaderParams(dirb, LongRange.NONE, SequenceMode.BIDIRECTIONAL);
    final ReaderParams b2 = DefaultReaderParamsTest.createDefaultReaderParams(dirb, LongRange.NONE, SequenceMode.UNIDIRECTIONAL);
    TestUtils.equalsHashTest(new ReaderParams[][] {{a1, a2}, {b1}, {b2}});
    a1.close();
    a2.close();
    b1.close();
    b2.close();
  }

  public void test0() throws IOException, ClassNotFoundException {
    final ReaderParams sp = getParams();
    final String dirStr = sp.directory().toString();
    assertTrue(dirStr.contains("test") && dirStr.contains("unit"));
    assertEquals(4, sp.maxLength());
    final SequencesReader reader = sp.reader();
    assertEquals(reader.path(), sp.directory());
    assertEquals(reader.type(), sp.mode().type());
    sp.close();
    assertTrue(sp.toString().startsWith("SequenceParams mode=BIDIRECTIONAL directory="));

  }

  public void testNotFound() throws IOException {
    final ByteArrayOutputStream bs = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bs));
    final ByteArrayOutputStream dls = new ByteArrayOutputStream();
    final PrintStream dlp = new PrintStream(dls);
    final DiagnosticListener listener = new DiagnosticListener() {
      @Override
      public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
        dlp.println(event.getMessage());
      }

      @Override
      public void close() {
      }
    };
    try {
      Diagnostic.addListener(listener);

      try {
        final ReaderParams sp = DefaultReaderParamsTest.createDefaultReaderParams(new File("nonexistantfile"), LongRange.NONE, SequenceMode.BIDIRECTIONAL);
        sp.reader();
        fail();
      } catch (final SlimException e) {
        e.logException();
        e.printErrorNoLog();
      }
      dlp.close();
      bs.close();

      final String actual = bs.toString();
      assertTrue(actual.contains("The specified SDF, \"nonexistantfile\", does not exist." + StringUtils.LS));
      dls.close();
      assertEquals("Error: The specified SDF, \"nonexistantfile\", does not exist." + StringUtils.LS, dls.toString());
    } finally {
      Diagnostic.removeListener(listener);
      Diagnostic.setLogStream();
    }
  }
}

