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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

/**
 */
public class SingletonValidatorTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SingletonCli();
  }

  public void test() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory("variancevalidator")) {
      final File sdf = new File(tempDir, "sdf");
      ReaderTestUtils.getDNADir(REF_SEQS, sdf);
      final File aFile = new File(tempDir, "anyFile");
      assertTrue(aFile.createNewFile());
      final File listFile = new File(tempDir, "list.listFile");
      FileUtils.stringToFile(aFile.getPath(), listFile);

      checkHandleFlagsErr("-o", tempDir.getPath(), "-t", "blah");
      checkHandleFlagsErr("-o", "blah", "-t", "blah");
      checkHandleFlagsErr("-o", "blah", "-t", sdf.getPath(), "--input-list-listFile", listFile.getPath());
      checkHandleFlagsErr("-o", "blah", "-t", sdf.getPath(), "blah", "-T", "0");
      checkHandleFlagsErr("-o", "blah", "-t", sdf.getPath(), "blah", "--ploidy", "diploid", "--sex", "male");
      checkHandleFlagsErr("-o", "blah", "-t", sdf.getPath(), "blah", "--ploidy", "diploid", "--pedigree", aFile.getPath());

      String err = checkHandleFlagsErr("-o", "blah", "-t", sdf.getPath(), aFile.getPath(), "--filter-ambiguity", "-1");
      assertTrue(err, err.contains("The specified flag \"--filter-ambiguity\" has invalid value \"-1\". It should be greater than or equal to \"0%\"."));

      err = checkHandleFlagsErr("-o", "blah", "-t", sdf.getPath(), aFile.getPath(), "--filter-ambiguity", "101%");
      assertTrue(err, err.contains("The specified flag \"--filter-ambiguity\" has invalid value \"101%\". It should be less than or equal to \"100%\"."));

      checkHandleFlagsOut("-o", "blah", "-t", sdf.getPath(), "-I", listFile.getPath());
    }
  }
}
