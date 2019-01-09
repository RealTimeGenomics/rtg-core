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
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SequenceParamsTest extends TestCase {

  private static final String SEQ_DNA_A = "" + ">x" + StringUtils.LS + "actg" + StringUtils.LS + ">y" + StringUtils.LS + "actg" + StringUtils.LS + ">z" + StringUtils.LS + "actg" + StringUtils.LS + ">a" + StringUtils.LS + "actg" + StringUtils.LS;

  private File mDir = null;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }
  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
    Diagnostic.setLogStream();
  }

  SequenceParams getParams() throws IOException {
    final File dir = getDir();
    return SequenceParams.builder().directory(dir).create();
  }

  SequenceParams getParams(final long start, final long end) throws IOException {
    final File dir = getDir();
    return SequenceParams.builder().directory(dir).mode(SequenceMode.BIDIRECTIONAL).readerRestriction(LongRange.NONE).useMemReader(false).loadNames(false).region(new HashingRegion(start, end)).create();
  }

  private File getDir() throws IOException {
    final File dir = FileUtils.createTempDir("SequenceParamsTest", "getDir", mDir);
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(SEQ_DNA_A.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams,
            new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, false);
    sequenceWriter.processSequences();
    return dir;
  }

  ISequenceParams getParams(final ReaderParams readerParams, final long start, final long end) {
    return new MockSequenceParams(readerParams, SequenceMode.BIDIRECTIONAL, start, end);
  }

  public void testSubSequence0() throws IOException {
    final ISequenceParams par = getParams(0, 4);
    final ISequenceParams sub = par.subSequence(new HashingRegion(1, 3));
    assertEquals(1, sub.region().getStart());
    assertEquals(3, sub.region().getEnd());
    assertFalse(par.reader() == sub.reader());
    par.close();
    sub.close();
  }

  public void testEquals() throws IOException {
    final File dira = getDir();
    final File dirb = getDir();
    final ISequenceParams a1 = SequenceParams.builder().directory(dira).mode(SequenceMode.BIDIRECTIONAL).region(new HashingRegion(0, 1)).create();
    final ISequenceParams a2 = SequenceParams.builder().directory(dira).mode(SequenceMode.BIDIRECTIONAL).region(new HashingRegion(0, 1)).create();
    final ISequenceParams b1 = SequenceParams.builder().directory(dirb).mode(SequenceMode.BIDIRECTIONAL).region(new HashingRegion(0, 1)).create();
    final ISequenceParams b2 = SequenceParams.builder().directory(dirb).mode(SequenceMode.UNIDIRECTIONAL).region(new HashingRegion(0, 1)).create();
    final ISequenceParams c = SequenceParams.builder().directory(dirb).mode(SequenceMode.BIDIRECTIONAL).region(new HashingRegion(1, 1)).create();
    final ISequenceParams d = SequenceParams.builder().directory(dirb).mode(SequenceMode.BIDIRECTIONAL).region(new HashingRegion(0, 2)).create();
    TestUtils.equalsHashTest(new ISequenceParams[][] {{a1, a2}, {b1}, {b2}, {c}, {d}});
    assertTrue(a1.toString().startsWith("SequenceParams mode=BIDIRECTIONAL region=[(0:-1), (1:-1)] directory="));
    a1.close();
    a2.close();
    b1.close();
    b2.close();
    c.close();
    d.close();
  }

  public void test0() throws IOException {
    final SequenceParams sp = getParams();
    assertEquals(LongRange.NONE, sp.readerRestriction());
    sp.integrity();
    final String dirStr = sp.directory().toString();
    assertTrue(dirStr.contains("test") && dirStr.contains("unit"));
    assertEquals(0, sp.region().getStart());
    assertEquals(4, sp.region().getEnd());
    assertEquals(4, sp.numberSequences());
    assertEquals(4, sp.maxLength());
    final SequencesReader reader = sp.reader();
    assertEquals(reader.path(), sp.directory());
    assertEquals(reader.type(), sp.mode().type());
    sp.close();
    assertTrue(sp.toString().startsWith("SequenceParams mode=BIDIRECTIONAL region=[(0:-1), (4:-1)] directory="));

  }

  public void test1() throws IOException {
    final SequenceParams sp = getParams(2, 3);
    sp.integrity();
    final String dirStr = sp.directory().toString();
    assertTrue(dirStr.contains("test") && dirStr.contains("unit"));
    assertEquals(2, sp.region().getStart());
    assertEquals(3, sp.region().getEnd());
    assertEquals(1, sp.numberSequences());
    final SequencesReader reader = sp.reader();
    assertEquals(reader.path(), sp.directory());
    assertEquals(reader.type(), sp.mode().type());
    sp.close();
    assertTrue(sp.toString(), sp.toString().startsWith("SequenceParams mode=BIDIRECTIONAL region=[(2:-1), (3:-1)] directory="));
  }

  public void test1a() throws IOException {
    final ISequenceParams sp0 = getParams();
    final ISequenceParams sp = getParams(sp0.readerParams(), 2, 3);
    final String dirStr = sp.directory().toString();
    assertTrue(dirStr.contains("test") && dirStr.contains("unit"));
    assertEquals(2, sp.region().getStart());
    assertEquals(3, sp.region().getEnd());
    assertEquals(1, sp.numberSequences());
    final SequencesReader reader = sp.reader();
    assertEquals(reader.path(), sp.directory());
    assertEquals(reader.type(), sp.mode().type());
    sp.close();
    assertTrue(sp.toString().startsWith("SequenceParams mode=BIDIRECTIONAL region=[(2:-1), (3:-1)] directory="));
  }

  public void test2() throws IOException {
    final ISequenceParams sp = getParams(2, 2);
    final String dirStr = sp.directory().toString();
    assertTrue(dirStr.contains("test") && dirStr.contains("unit"));
    assertEquals(2, sp.region().getStart());
    assertEquals(2, sp.region().getEnd());
    assertEquals(0, sp.numberSequences());
    final SequencesReader reader = sp.reader();
    assertEquals(reader.path(), sp.directory());
    assertEquals(reader.type(), sp.mode().type());
    sp.close();
    assertTrue(sp.toString().startsWith("SequenceParams mode=BIDIRECTIONAL region=[(2:-1), (2:-1)] directory="));
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
        final ISequenceParams sp = SequenceParams.builder().directory(new File("nonexistantfile")).create();
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

