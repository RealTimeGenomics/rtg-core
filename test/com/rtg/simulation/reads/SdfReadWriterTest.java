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

package com.rtg.simulation.reads;

import java.io.File;

import com.rtg.reader.PrereadType;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.AbstractTempFileHandler;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class SdfReadWriterTest extends AbstractTempFileHandler {

  public void testPaired() throws Exception {
    Diagnostic.setLogStream();
    final File sdf = new File(mTempDir, "sdf");
    SdfReadWriter w = new SdfReadWriter(sdf, true, PrereadType.SOLEXA, true, true);
    assertEquals(0, w.readsWritten());
    try {
      w.writeRead("read", new byte[] {1, 2, 3, 4}, new byte[] {1, 2, 3, 4}, 4);
      fail();
    } catch (final IllegalStateException e) {
      // ok
    }
    w.writeLeftRead("read", new byte[] {1, 2, 3, 4}, new byte[] {1, 2, 3, 4}, 4);
    w.writeRightRead("read", new byte[] {1, 2, 3, 4}, new byte[] {1, 2, 3, 4}, 4);
    w.close();
    assertEquals(1, w.readsWritten());
    final File fasta = new File(mTempDir, "f.fasta.gz");
    assertEquals(0, new Sdf2Fasta().mainInit(new String[] {"-i", sdf.getPath(), "-o", fasta.getPath()}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
    assertEquals(">0 read" + StringUtils.LS + "ACGT" + StringUtils.LS, FileHelper.gzFileToString(new File(mTempDir, "f_1.fasta.gz")));
    assertEquals(">0 read" + StringUtils.LS + "ACGT" + StringUtils.LS, FileHelper.gzFileToString(new File(mTempDir, "f_2.fasta.gz")));
  }

  public void testSingle() throws Exception {
    Diagnostic.setLogStream();
    final File sdf = new File(mTempDir, "sdf");
    SdfReadWriter w = new SdfReadWriter(sdf, false, PrereadType.SOLEXA, true, true);
    try {
      w.writeLeftRead("read", new byte[] {1, 2, 3, 4}, new byte[] {1, 2, 3, 4}, 4);
      fail();
    } catch (final IllegalStateException e) {
      // ok
    }
    try {
      w.writeRightRead("read", new byte[] {1, 2, 3, 4}, new byte[] {1, 2, 3, 4}, 4);
      fail();
    } catch (final IllegalStateException e) {
      // ok
    }
    w.writeRead("read", new byte[] {1, 2, 3, 4}, new byte[] {1, 2, 3, 4}, 4);
    w.close();
    assertEquals(1, w.readsWritten());
    final File fasta = new File(mTempDir, "f.fasta.gz");
    assertEquals(0, new Sdf2Fasta().mainInit(new String[] {"-i", sdf.getPath(), "-o", fasta.getPath()}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
    assertEquals(">0 read" + StringUtils.LS + "ACGT" + StringUtils.LS, FileHelper.gzFileToString(new File(mTempDir, "f.fasta.gz")));
  }
}
