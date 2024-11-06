/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
