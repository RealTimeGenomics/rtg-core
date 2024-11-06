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
