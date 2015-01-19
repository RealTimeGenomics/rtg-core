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

package com.rtg.variant.bayes.multisample.singleton;

import static com.rtg.sam.SharedSamConstants.REF_SEQS;

import java.io.File;
import java.io.IOException;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class SingletonValidatorTest extends TestCase {

  public void test() throws IOException {
    final MemoryPrintStream log = new MemoryPrintStream();
    Diagnostic.setLogStream(log.printStream());
    try (final TestDirectory tempDir = new TestDirectory("variancevalidator")) {
      final File sdf = new File(tempDir, "sdf");
      ReaderTestUtils.getDNADir(REF_SEQS, sdf);
      final File aFile = new File(tempDir, "anyFile");
      assertTrue(aFile.createNewFile());
      final File listFile = new File(tempDir, "list.listFile");
      FileUtils.stringToFile(aFile.getPath(), listFile);
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final CFlags flags = new CFlags("", out.printStream(), err.printStream());
      new SingletonCli().initLocalFlags(flags);

      String[] inputFlags = {
          "-o", tempDir.getPath(),
          "-t", "blah",
      };
      assertFalse(flags.setFlags(inputFlags));

      inputFlags = new String[] {
          "-o", "blah",
          "-t", "blah",
      };
      assertFalse(flags.setFlags(inputFlags));

      inputFlags = new String[] {
          "-o", "blah",
          "-t", sdf.getPath(),
          "--input-list-listFile", listFile.getPath()
      };
      assertFalse(flags.setFlags(inputFlags));

      inputFlags = new String[] {
          "-o", "blah",
          "-t", sdf.getPath(),
          "blah",
          "-T", "0"
      };
      assertFalse(flags.setFlags(inputFlags));

      inputFlags = new String[] {
          "-o", "blah",
          "-t", sdf.getPath(),
          aFile.getPath(),
          "--filter-ambiguity", "-1"
      };
      assertFalse(flags.setFlags(inputFlags));
      assertTrue(log.toString(), log.toString().contains("The specified flag \"--filter-ambiguity\" has invalid value \"-1\". It should be greater than or equal to \"0%\"."));

      inputFlags = new String[] {
          "-o", "blah",
          "-t", sdf.getPath(),
          aFile.getPath(),
          "--filter-ambiguity", "101%"
      };
      assertFalse(flags.setFlags(inputFlags));
      assertTrue(log.toString(), log.toString().contains("The specified flag \"--filter-ambiguity\" has invalid value \"101%\". It should be less than or equal to \"100%\"."));

      inputFlags = new String[] {
          "-o", "blah",
          "-t", sdf.getPath(),
          "-I", listFile.getPath()
      };
      if (!flags.setFlags(inputFlags)) {
        fail(flags.getParseMessage());
      }

      //      System.err.println("log:");
      //      System.err.println(log.toString());
      //      System.err.println("err:");
      //      System.err.println(err.toString());
      //      System.err.println("out:");
      //      System.err.println(out.toString());
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
