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

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorEvent;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Tests for corresponding class.
 */
public class Sdf2FastqTest extends AbstractCliTest {


  @Override
  protected AbstractCli getCli() {
    return new Sdf2Fastq();
  }

  public void testHelp() {
    checkHelp("output filename (extension added if not present)",
      "SDF containing sequences"
    );
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
    try (TestDirectory dir = new TestDirectory()) {
      final File sdf = ReaderTestUtils.getDNAFastqDir(FULL_NAME_DATA, new File(dir, "sdf"), false);
      final File fasta = new File(dir, "fs.fastq.gz");
      checkMainInitOk("-i", sdf.getPath(), "-o", fasta.getPath());
      assertEquals(FULL_NAME_DATA, FileHelper.gzFileToString(fasta));
    }
  }
}
