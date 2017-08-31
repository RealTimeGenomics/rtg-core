/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.sv;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class FusionFilterTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new FusionFilter();
  }

  public void testErrorMessages() {
    checkHelp("Filters a VCF containing structural variant", "acceptor gene regions", "donor gene regions", "fusion-scores");
    TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-a", "acceptor.bed", "-d", "donor.bed", "-o", "-", "-i", "in.vcf"), "does not exist");
  }

  public void testSimple() throws IOException, UnindexableDataException {
    endToEnd("fusion", "fusion");
  }

  protected void endToEnd(String harnessId, String resultsId, String... args) throws IOException, UnindexableDataException {
    try (TestDirectory dir = new TestDirectory("fusionfilter-nano")) {
      final File inVcf = new File(dir, "input.vcf.gz");
      FileHelper.stringToGzFile(mNano.loadReference(harnessId + "_in.vcf"), inVcf);
      new TabixIndexer(inVcf).saveVcfIndex();
      final File donorBed = new File(dir, "donor.bed");
      FileUtils.stringToFile(mNano.loadReference(harnessId + "_in_donor.bed"), donorBed);
      final File acceptorBed = new File(dir, "acceptor.bed");
      FileUtils.stringToFile(mNano.loadReference(harnessId + "_in_acceptor.bed"), acceptorBed);
      final File scoreTsv = new File(dir, "fusion-scores.tsv");
      FileUtils.stringToFile(mNano.loadReference(harnessId + "_in_score.tsv"), scoreTsv);

      final File outVcf = new File(dir, "filtered.vcf");

      final String[] fullArgs = Utils.append(args, "-a", acceptorBed.getPath(), "-d", donorBed.getPath(), "-s", scoreTsv.getPath(), "-i", inVcf.getPath(), "-o", outVcf.getPath(), "-Z");
      final MainResult res = MainResult.run(getCli(), fullArgs);
      assertEquals(res.err(), 0, res.rc());

      final String content = FileUtils.fileToString(outVcf);
      mNano.check(resultsId + "_out_" + outVcf.getName(), TestUtils.sanitizeVcfHeader(content));
    }
  }
}
