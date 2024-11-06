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

package com.rtg.variant.cnv.segment;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Tests the corresponding class.
 */
public class CnvPonBuildCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CnvPonBuildCli();
  }

  public void testHelp() {
    checkHelp("Construct a normalized coverage sample",
      "BED output file",
      "reference genome",
      "coverage BED file. Must be specified 1 or more times",
      "do not gzip the output",
      "Utility",
      "File Input/Output"
      );
  }

  public void testErrors() {
    TestUtils.containsAll(checkHandleFlagsErr(), "Error: You must provide values for -o FILE -t SDF FILE+");
  }

  public void testValid() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File coverage = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/coverage.bed", new File(dir, "coverage.bed"));
      final File output = new File(dir, "panel.bed");
      // Slightly G rich reference so G+C normalization actually does something
      final String ref = ">1\n" + StringUtils.repeat("ACGGT", 25000 / "ACGGT".length());
      final File sdf = new File(dir, "sdf");
      ReaderTestUtils.getDNADir(ref, sdf);
      final MainResult result = checkMainInit("-Z", "-o", output.getPath(), "-t", sdf.getPath(), "--label-column-name=label", coverage.getPath());
      final String actual = StringUtils.grepMinusV(FileHelper.fileToString(output), "^#");
      mNano.check("expected.panel.bed", actual);
      TestUtils.containsAll(result.out(), "Normalizing and G+C correcting", "coverage.bed");
    }
  }

}
