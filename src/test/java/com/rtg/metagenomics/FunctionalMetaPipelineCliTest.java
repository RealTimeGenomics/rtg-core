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

package com.rtg.metagenomics;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.io.TestDirectory;

/**
 */
public class FunctionalMetaPipelineCliTest extends AbstractCliTest {


  public void testFlags() throws IOException {
    try (TestDirectory tmp = new TestDirectory()) {
      final File protein = new File(tmp, "protein");
      final File dna = new File(tmp, "dna");
      final File readSdf = new File(tmp, "dna");
      ReaderTestUtils.getReaderProtein(">a" + LS + "GGATASDCASSZZXCVB" + LS, protein);
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGTTTAGA" + LS, dna, new SdfId());
      checkHandleFlagsOut("--input", readSdf.getPath(), "--filter", dna.getPath(), "--protein", protein.getPath(), "--output", new File(tmp, "output").getPath());
    }
  }

  @Override
  protected AbstractCli getCli() {
    return new FunctionalMetaPipelineCli();
  }
}
