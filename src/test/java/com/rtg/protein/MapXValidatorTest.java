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
package com.rtg.protein;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Test for corresponding class
 */
public class MapXValidatorTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new MapXCli();
  }

  private static final String TEMPLATE_DNA = "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGAAATTATAACCACGACGCAGCAGACGCAG";
  private static final String TEMPLATE_PROTEIN = TestUtils.dnaToProtein(TEMPLATE_DNA);
  private static final String TEMPLATE_FASTA = ">templateName" + LS + TEMPLATE_PROTEIN + LS;

  private static final String[] READS_ONE_INDEL = {
    "TGGCG" + "ACG" + "CAAAAACAGAAAGTCGAAAAAAAATCAA",
    "AATGGCGCAAAAACAGAAAGTCATGGAAAAAAAATC",
    "ATGGCGCAAAAACACGCGAAAGTCGAAAAAAAATCA",
    "TGGCGCAAAAACAGAAAGTCGAAATTTAAAAATCAA",
    DnaUtils.reverseComplement("AAATGGCGCAAAAAACCCAGAAAGTCGAAAAAAAAT"),
    DnaUtils.reverseComplement("AATGGCGCAAAAACAGAAAGTTGACGAAAAAAAATC"),
    DnaUtils.reverseComplement("ATGGCGCAAAAACAGAAAGTCGGGGAAAAAAAATCA"),
    DnaUtils.reverseComplement("TGGCGCAAAAACAGAAAGTCGACACAAAATCACGAA"),
  "AAAAAAATCAAAGAAATTATAACCACGACAAAGCAG"};

  private static final String READS_FASTA_ONE_INDEL;
  static {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < READS_ONE_INDEL.length; ++i) {
      sb.append(">testRead").append(i).append(LS).append(READS_ONE_INDEL[i]).append(LS);
    }
    READS_FASTA_ONE_INDEL = sb.toString();
  }

  private void assertStringContains(String base, String sub) {
    assertTrue(base, base.contains(sub));
  }

  public void testMapXValidator() throws IOException {
    try (final TestDirectory dir = new TestDirectory("mapx")) {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      final File left = new File(reads, "left");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(READS_FASTA_ONE_INDEL, left, null).close();

      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=noblosum"), "Error: Invalid value \"noblosum\" for flag --matrix=noblosum.");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-l", left.getPath(), "-r", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1"), "Unknown flag -l");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "-n", "251"), "--max-top-results must be in the range [1, 250]");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "-n", "0"), "--max-top-results must be in the range [1, 250]");
      assertStringContains(checkHandleFlagsErr("-t", "blaht", "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "13"), "The specified SDF, \"blaht\", does not exist");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-identity", "-1"), " must be in the range [0, 100]");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-identity", "101"), "must be in the range [0, 100]");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "101", "--min-bit-score", "5"), "Cannot set both");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "-0.01", "--min-bit-score", "-1"), "Cannot set both");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "-0.01"), "at least 0.0");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-bit-score", "-0.01"), "at least 0.0");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-bit-score", "0.5", "-n", "0"), "--max-top-results must be in the range [1, 250]");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "0.5", "-n", "0"), "--max-top-results must be in the range [1, 250]");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "13"), "--word must be in the range [1, 12]");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "-1", "-b", "0", "-w", "12"), "--mismatches must be at least 0");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "0", "-b", "-1", "-w", "12"), "--gaps must be at least 0");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-c", "0"), "--gap-length must be at least 1");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=blosum45", "--step", "-1"), "Unknown flag --step");

      checkMainInitOk("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=blosum45");
      checkMainInitOk("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", new File(dir, "output2").getPath(), "-w", "4", "-T", "1", "--matrix=blosum80");
      FileHelper.deleteAll(output);
      checkMainInitOk("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=blosum45", "-f", "topn");

      final File right = new File(reads, "right");
      assertTrue(right.mkdir());
      ReaderTestUtils.getReaderDNA(READS_FASTA_ONE_INDEL, right, null).close();

      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", new File(dir, "output3").toString(), "-w", "4", "-T", "1", "--matrix=blosum45"), "Paired end data not supported");
    }
  }

}
