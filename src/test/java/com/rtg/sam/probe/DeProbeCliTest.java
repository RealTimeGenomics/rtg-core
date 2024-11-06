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
