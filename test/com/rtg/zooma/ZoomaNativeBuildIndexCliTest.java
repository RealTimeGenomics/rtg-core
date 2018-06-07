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
package com.rtg.zooma;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.SpawnJvm;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.RandomDna;

/**
 */
public class ZoomaNativeBuildIndexCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ZoomaNativeBuildIndexCli();
  }

  // called by main() in the spawned JVM
  private void runActualTest() throws IOException {
    setUp();
    try (final TestDirectory dir = new TestDirectory()) {
      final String reference = RandomDna.random(1000);
      final File referenceFile = new File(dir, "reference.fasta");
      final File indexFile = new File(dir, "reference.bin");
      FileUtils.stringToFile(">reference\n" + reference + "\n", referenceFile);
      checkMainInitOk("-i", referenceFile.toString(), "-o", indexFile.toString());
    }
  }

  // invoked by spawn in test() method
  public static void main(final String[] args) throws IOException {
    new ZoomaNativeBuildIndexCliTest().runActualTest();
  }

  public void test() throws IOException {
    if (NativeZooma.isEnabled()) {
      final Process p = SpawnJvm.spawn(ZoomaNativeBuildIndexCliTest.class.getName());
      final String out = FileUtils.streamToString(p.getErrorStream());
      TestUtils.containsAll(out, "max nucleotides to hash: 1000", "saveindex time");
      assertEquals("", FileUtils.streamToString(p.getInputStream()));
      p.getOutputStream().close();
      p.getInputStream().close();
      p.getErrorStream().close();
    } else {
      System.err.println("skipping zooma test (no native lib)");
    }
  }
}
