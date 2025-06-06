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
package com.rtg.variant.sv;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SvToolCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SvToolCli();
  }

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -t SDF";

  public void testErrorFlags() {
    checkParamsError(new String[] {}, EXP_F1);
  }

  private void checkParamsError(final String[] args, final String... exp) {
    TestUtils.containsAllUnwrapped(checkHandleFlagsErr(args), exp);
  }

  private static final String EXP_P1 = "The specified SDF, \"blah\", does not exist";
  private static final String[] EXP_P3 = {"The specified file, ", ", is not an SDF"};

  public void testErrorParams() throws IOException, InvalidParamsException {
    try (final TestDirectory tempDir = new TestDirectory()) {
      final File out = new File(tempDir, "out");
      checkParamsError(new String[]{"-o", out.getPath(), "-t", "blah", "blah.alignments"}, EXP_P1);
      FileHelper.deleteAll(out);
      final File tmp = File.createTempFile("test", "tmp", tempDir);
      final File aln = new File(tempDir, "blah.alignments");
      checkParamsError(new String[]{"-o", out.getPath(), "-t", tmp.getPath(), aln.getPath()}, EXP_P3);
      final File gen = FileUtils.createTempDir("test", "genome", tempDir);
      assertTrue(aln.createNewFile());
      final File rg = new File(tempDir, "rgfile.txt");
      assertTrue(rg.createNewFile());
      checkHandleFlagsOut("-o", out.getPath(), "-t", gen.getPath(), aln.getPath(), "-r", rg.getPath());

      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-b", "0"), "Expected a positive integer for parameter \"Xbin-size\"");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-s", "0"), "Expected a positive integer for parameter \"step\"");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-f", "0"), "Expected a positive integer for parameter \"fine-step\"");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", out.getPath(), "-t", gen.getPath(), "-r", rg.getPath(), aln.getPath(), "-s", "1"), "Parameter \"fine-step\" should be smaller than or equal to parameter \"step\"");
    }
  }

  public void testInitParams() {
    checkHelp("rtg sv"
        , "Analyses SAM records to determine the location of structural variants."
        , "[OPTION]... -o DIR -t SDF -I FILE -R FILE"
        , "[OPTION]... -o DIR -t SDF -r FILE FILE+"
        , "simple-signals", "if set, also output simple signals"
        , "-f,", "--fine-step=INT", "step size in interesting regions (Default is 10)"
        , "-s,", "--step=INT", "step size (Default is 100)"
        );
    checkExtendedHelp("rtg sv"
        , "-b,", "--Xbin-size=INT", "bin size used by simple signals (Default is 10)"
        , "--Xcorrections=FILE", "file containing per position corrections"
        , "--Xheterozygous", "if set, also include heterozygous bayesian models"
        );
  }
}
