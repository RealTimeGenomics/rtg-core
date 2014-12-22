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
 * Tests for corresponding class.
 */
public class Sdf2FastqTest extends TestCase {

  public Sdf2FastqTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(Sdf2FastqTest.class);
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void testHelp() {
    final CFlags flags = new CFlags();
    new Sdf2Fastq().initFlags(flags);
    TestCFlags.check(flags,
                     "output filename (extension added if not present)",
                     "SDF containing sequences"
                     );
  }

  public void testProcess() throws IOException, InvalidParamsException {
    ByteArrayOutputStream b = new ByteArrayOutputStream();
    LineWriter p = new LineWriter(new OutputStreamWriter(b));
    Sdf2Fastq.process(new ArraySequencesReader(new byte[][] {{1}}, new byte[][] {{0}}), p, false, 0, -1);
    p.flush();
    assertEquals("@sequence 0" + StringUtils.LS + "A" + StringUtils.LS + "+sequence 0" + StringUtils.LS + "!" + StringUtils.LS,
        b.toString());

    b = new ByteArrayOutputStream();
    p = new LineWriter(new OutputStreamWriter(b));
    Sdf2Fastq.process(new ArraySequencesReader(new byte[][] {{1, 2, 3, 4}}, new byte[][] {{0, 42, 85, 86}}), p, false, 2, -1);
    p.flush();
    assertEquals("@sequence 0" + StringUtils.LS + "AC" + StringUtils.LS + "GT" + StringUtils.LS + "+sequence 0" + StringUtils.LS + "!K" + StringUtils.LS + "vw" + StringUtils.LS,
        b.toString());
  }

  public void testValidator() {
    final int[] blah = new int[1];
    final DiagnosticListener dl = new DiagnosticListener() {
      @Override
      public void handleDiagnosticEvent(DiagnosticEvent<?> event) {
        if (event instanceof ErrorEvent) {
          assertEquals("Error: Expected a nonnegative integer for parameter \"line-length\".", event.getMessage());
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
      final Sdf2Fastq ptfq = new Sdf2Fastq();
      assertEquals(1, ptfq.mainInit(new String[] {"-i", "blah", "-o", "blaho", "-l", "-3"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      assertEquals(1, blah[0]);
    } finally {
      Diagnostic.removeListener(dl);
    }
  }

  static final String FULL_NAME_DATA = ""
          + "@name suffix" + StringUtils.LS
          + "ACGTCG" + StringUtils.LS
          + "+name suffix" + StringUtils.LS
          + "123456" + StringUtils.LS
          + "@second suffixes" + StringUtils.LS
          + "ACGGGT" + StringUtils.LS
          + "+second suffixes" + StringUtils.LS
          + "123456" + StringUtils.LS;

  public void testFullName() throws IOException {
    final File dir = FileUtils.createTempDir("testsdf2fasta", "fullname");
    try {
      final File sdf = ReaderTestUtils.getDNAFastqDir(FULL_NAME_DATA, new File(dir, "sdf"), false);
      final File fasta = new File(dir, "fs.fastq.gz");
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final int code = new Sdf2Fastq().mainInit(new String[] {"-i", sdf.getPath(), "-o", fasta.getPath()}, out.outputStream(), err.printStream());
      assertEquals(err.toString(), 0, code);
      final String outStr = FileHelper.gzFileToString(fasta);
      assertEquals(FULL_NAME_DATA, outStr);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
