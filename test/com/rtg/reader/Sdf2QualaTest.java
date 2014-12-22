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
package com.rtg.reader;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorEvent;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for the corresponding class.
 */
public class Sdf2QualaTest extends TestCase {

  public Sdf2QualaTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(Sdf2QualaTest.class);
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
  }

  public void testProcess1() throws IOException, InvalidParamsException {
    final ByteArrayOutputStream seqStream = new ByteArrayOutputStream();
    final LineWriter seqWriter = new LineWriter(new OutputStreamWriter(seqStream));
    final ByteArrayOutputStream qualStream = new ByteArrayOutputStream();
    final LineWriter qualWriter = new LineWriter(new OutputStreamWriter(qualStream));
    Sdf2Quala.process(new ArraySequencesReader(new byte[][] {{1}}, new byte[][] {{0}}), seqWriter, qualWriter, false, (byte) -1);
    seqWriter.close();
    qualWriter.close();
    assertEquals(">sequence 0" + StringUtils.LS + "A" + StringUtils.LS, seqStream.toString());
    assertEquals(">sequence 0" + StringUtils.LS + "0" + StringUtils.LS, qualStream.toString());
  }

  public void testProcess2() throws IOException, InvalidParamsException {
    final ByteArrayOutputStream seqStream = new ByteArrayOutputStream();
    final LineWriter seqWriter = new LineWriter(new OutputStreamWriter(seqStream));
    final ByteArrayOutputStream qualStream = new ByteArrayOutputStream();
    final LineWriter qualWriter = new LineWriter(new OutputStreamWriter(qualStream));
    Sdf2Quala.process(new ArraySequencesReader(new byte[][] {{1, 2, 3, 4}}, new byte[][] {{0, 42, 85, 86}}), seqWriter, qualWriter, false, (byte) -1);
    seqWriter.close();
    qualWriter.close();
    assertEquals(">sequence 0" + StringUtils.LS + "ACGT" + StringUtils.LS, seqStream.toString());
    assertEquals(">sequence 0" + StringUtils.LS + "0 42 85 86" + StringUtils.LS, qualStream.toString());
  }

  private static final String FULL_NAME_SEQ_DATA = ""
          + ">name suffix" + StringUtils.LS
          + "ACGTCG" + StringUtils.LS
          + ">second suffixes" + StringUtils.LS
          + "ACGGGT" + StringUtils.LS
          ;

  private static final String FULL_NAME_QUAL_DATA = ""
          + ">name suffix" + StringUtils.LS
          + "16 17 18 19 20 21" + StringUtils.LS
          + ">second suffixes" + StringUtils.LS
          + "16 17 18 19 20 21" + StringUtils.LS
          ;

  public void testFullName() throws IOException {
    final File dir = FileUtils.createTempDir("testsdf2fasta", "fullname");
    try {
      final File sdf = ReaderTestUtils.getDNAFastqDir(Sdf2FastqTest.FULL_NAME_DATA, new File(dir, "sdf"), false);
      final File base = new File(dir, "fs");
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final int code = new Sdf2Quala().mainInit(new String[] {"-i", sdf.getPath(), "-o", base.getPath()}, out.outputStream(), err.printStream());
      assertEquals(err.toString(), 0, code);
      assertEquals(FULL_NAME_SEQ_DATA, FileHelper.gzFileToString(new File(base + ".fasta.gz")));
      assertEquals(FULL_NAME_QUAL_DATA, FileHelper.gzFileToString(new File(base + ".quala.gz")));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testHelp() {
    final CFlags flags = new CFlags();
    Sdf2Quala.initFlags(flags);
    TestCFlags.check(flags,
                     "basename for output files (extensions will be added)",
                     "SDF containing sequences"
                     );
  }

  private void checkQualityBadValue(final int qual) {
    final int[] blah = new int[1];
    final DiagnosticListener dl = new DiagnosticListener() {
      @Override
      public void handleDiagnosticEvent(DiagnosticEvent<?> event) {
        if (event instanceof ErrorEvent) {
          final String msg = event.getMessage();
          assertTrue(msg.startsWith("Error: The specified flag \"default-quality\" has invalid value \"" + qual + "\". It should be"));
          blah[0] += 1;
        } else {
          fail();
        }
      }
      @Override
      public void close() {
      }
    };
    Diagnostic.addListener(dl);
    try {
      final Sdf2Quala ptfq = new Sdf2Quala();
      assertEquals(1, ptfq.mainInit(new String[] {"-i", "blah", "-o", "blaho", "-q", String.valueOf(qual)}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      assertEquals(1, blah[0]);
    } finally {
      Diagnostic.removeListener(dl);
    }
  }

  public void testValidator() {
    checkQualityBadValue(-1);
    checkQualityBadValue(64);
  }
}
