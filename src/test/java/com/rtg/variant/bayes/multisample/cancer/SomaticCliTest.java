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
package com.rtg.variant.bayes.multisample.cancer;

import static com.rtg.util.StringUtils.LS;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

/**
 */
public class SomaticCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SomaticCli();
  }

  private static final String EXP_F1 = "Error: You must provide values for -o DIR -t SDF" + LS;

  public void testErrorF1() {
    assertTrue(checkHandleFlagsErr().contains(EXP_F1));
  }

  public void testInitParams() {
    checkHelp("somatic [OPTION]... -o DIR -t SDF --contamination FLOAT --derived STRING --original STRING FILE+",
      "sample identifier used in read groups for original",
      "threshold for ambiguity",
      "include gain of reference somatic calls in output VCF",
      "estimated fraction of contamination"
    );
  }

  public void testValidator() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File tmpFile = FileUtils.stringToFile("original-derived TEST cancer contamination=0.13", new File(tmpDir, "relations.relations"));
      final File in = new File(tmpDir, "alignments.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM_CANCER, in);
      new TabixIndexer(in, new File(tmpDir, "alignments.sam.gz.tbi")).saveSamIndex();

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getDNADir(">g1\nacgtacgtacgtacgtacgt", template);
      final File outDir = new File(tmpDir, "out");
      String err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", outDir.getPath());
      TestUtils.containsAllUnwrapped(err, "Error: The specified SDF, \"" + outDir.getPath() + "\", does not exist.");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--original", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "All of --derived, --original, and --contamination must be specified");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--derived", "blah", in.getPath());
      TestUtils.containsAllUnwrapped(err, "All of --derived, --original, and --contamination must be specified");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "All of --derived, --original, and --contamination must be specified");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--original", "blah", in.getPath(), "-r", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --Xpedigree in conjunction with --derived, --original, or --contamination");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--derived", "blah", in.getPath(), "-r", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --Xpedigree in conjunction with --derived, --original, or --contamination");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--contamination", "0.5", in.getPath(), "-r", tmpFile.getPath());
      TestUtils.containsAllUnwrapped(err, "Cannot use --Xpedigree in conjunction with --derived, --original, or --contamination");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--sex", "male", "--original", "blah", "--derived", "blah", "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Original and derived must be different samples");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--sex", "male", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "Sex-specific processing was specified but " + template.getPath() + " is missing a 'reference.txt'");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--somatic", "1", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "--somatic must be in the range (0.0, 1.0)");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--somatic", "0", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "--somatic must be in the range (0.0, 1.0)");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "-0.5", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "--loh must be in the range [0.0, 1.0]");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "1.5", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      TestUtils.containsAllUnwrapped(err, "--loh must be in the range [0.0, 1.0]");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--contamination", "-0.5", "--original", "foo", "--derived", "bar", in.getPath());
      TestUtils.containsAllUnwrapped(err, "--contamination must be in the range [0.0, 1.0)");

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--contamination", "1.5", "--original", "foo", "--derived", "bar", in.getPath());
      TestUtils.containsAllUnwrapped(err, "--contamination must be in the range [0.0, 1.0)");


      checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "0", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "1", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());

      checkMainInitOk("-o", outDir.getPath(), "-t", template.getPath(), "--Xpedigree", tmpFile.getPath(), in.getPath(), in.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
    }
  }

  public void testPruneHypotheseFlag() {
    final CFlags flags = new CFlags();
    final SomaticCli cli = new SomaticCli();
    cli.initLocalFlags(flags);
    //This flag has to be true to prevent hypothesis overrun
    assertTrue((Boolean) flags.getValue("Xprune-hypotheses"));
  }
}
