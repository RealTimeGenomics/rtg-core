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
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.MainResult;
import com.rtg.sam.SamUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class DeProbeCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new DeProbeCli();
  }

  public void testValidator() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File probes = FileHelper.resourceToFile("com/rtg/sam/probe/resources/probes.bed", new File(dir, "probes.bed"));
      final File alignments = FileHelper.resourceToFile("com/rtg/sam/probe/resources/alignments.sam", new File(dir, "alignments.sam"));
      final File output = new File(dir, "output");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-o", output.getPath(), "-b", probes.getPath(), "--tolerance", "3"), "input files");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr(alignments.getPath(), "-o", output.getPath(), "-b", output.getPath(), "--tolerance", "3"), "--probe-bed", "does not exist");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr(alignments.getPath(), "-o", output.getPath(), "-b", probes.getPath(), "--tolerance", "-3"), "--tolerance", "at least 0");
      checkHandleFlagsOut(alignments.getPath(), "-o", output.getPath(), "-b", probes.getPath(), "--tolerance", "3");
    }
  }

  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File probes = FileHelper.resourceToFile("com/rtg/sam/probe/resources/probes.bed", new File(dir, "probes.bed"));
      final File alignments = FileHelper.resourceToFile("com/rtg/sam/probe/resources/alignments.sam", new File(dir, "alignments.sam"));
      final File output = new File(dir, "output");
      final File bamFile = new File(output, DeProbeCli.ALIGNMENT_FILE_NAME);
      final MainResult result = checkMainInit(alignments.getPath(), "-o", output.getPath(), "-b", probes.getPath(), "--tolerance", "3");
      mNano.check("expected.stripped.sam", SamUtils.bamToString(bamFile));
      mNano.check("expected.probes.tsv", FileUtils.fileToString(new File(output, DeProbeCli.PROBE_OFFSET_TABLE_FILE)));
      mNano.check("expected.cigar_ops.tsv", FileUtils.fileToString(new File(output, DeProbeCli.CIGAR_OP_TABLE_FILE)));
      mNano.check("expected.on_target.tsv", FileUtils.fileToString(new File(output, DeProbeCli.ON_TARGET_SUMMARY_FILE)));
      mNano.check("expected.out.txt", FileUtils.fileToString(new File(output, CommonFlags.SUMMARY_FILE)));
      mNano.check("expected.out.txt", result.out());
      mNano.check("expected.summary.tsv", FileUtils.fileToString(new File(output, DeProbeCli.PROBE_SUMMARY_FILE)));
      mNano.check("positive_strand_probe_counts.bed", FileHelper.zipFileToString(new File(output, DeProbeCli.POS_COUNTS_NAME)));
      mNano.check("negative_strand_probe_counts.bed", FileHelper.zipFileToString(new File(output, DeProbeCli.NEG_COUNTS_NAME)));
    }
  }

}
