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

package com.rtg.sam.probe;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class BamStripProbesTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new BamStripProbes();
  }

  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File probes = FileHelper.resourceToFile("com/rtg/sam/probe/resources/probes.bed", new File(dir, "probes.bed"));
      final File alignments = FileHelper.resourceToFile("com/rtg/sam/probe/resources/alignments.sam", new File(dir, "alignments.sam"));
      final File output = new File(dir, "output.sam");
      final MainResult result = checkMainInit("-i", alignments.getPath(), "-o", output.getPath(), "-b", probes.getPath(), "-Z");
      mNano.check("expected.stripped.sam", FileHelper.fileToString(output));
      mNano.check("expected.out.text", result.out());
    }
  }

}
