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
package com.rtg.variant.bayes.multisample.cancer;

import static com.rtg.util.StringUtils.LS;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

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
      assertTrue(err, err.contains("Error: The specified SDF, \"" + outDir.getPath() + "\", does not exist."));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--original", "blah", in.getPath());
      assertTrue(err, err.contains("All of --derived, --original, and --contamination must be specified"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--derived", "blah", in.getPath());
      assertTrue(err, err.contains("All of --derived, --original, and --contamination must be specified"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("All of --derived, --original, and --contamination must be specified"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--original", "blah", in.getPath(), "-r", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --Xpedigree in conjunction with --derived, --original, or --contamination"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--derived", "blah", in.getPath(), "-r", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --Xpedigree in conjunction with --derived, --original, or --contamination"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--contamination", "0.5", in.getPath(), "-r", tmpFile.getPath());
      assertTrue(err, err.contains("Cannot use --Xpedigree in conjunction with --derived, --original, or --contamination"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--sex", "male", "--original", "blah", "--derived", "blah", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("Original and derived must be different samples"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--sex", "male", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("Sex-specific processing was specified but " + template.getPath() + " is missing a 'reference.txt'"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--somatic", "1", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("--somatic should be a probability 0<s<1"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--somatic", "0", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("--somatic should be a probability 0<s<1"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "-0.5", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("--loh should be a probability 0<=p<=1"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "1.5", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("--loh should be a probability 0<=p<=1"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--Xcontrary-probability", "-0.5", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("--Xcontrary-probability should be a probability 0<p<=1"));

      err = checkHandleFlagsErr("-o", outDir.getPath(), "-t", template.getPath(), "--Xcontrary-probability", "1.5", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      assertTrue(err, err.contains("--Xcontrary-probability should be a probability 0<p<=1"));

      checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "0", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());
      checkHandleFlagsOut("-o", outDir.getPath(), "-t", template.getPath(), "--loh", "1", "--original", "foo", "--derived", "bar", "--contamination", "0.5", in.getPath());

      checkMainInitOk("-o", outDir.getPath(), "-t", template.getPath(), "--Xpedigree", tmpFile.getPath(), in.getPath(), in.getPath());
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
