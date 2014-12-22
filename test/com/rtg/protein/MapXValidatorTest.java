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
    for (int i = 0; i < READS_ONE_INDEL.length; i++) {
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

      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=noblosum"), "Error: Invalid value \"noblosum\" for \"--matrix=noblosum\".");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-l", left.getPath(), "-r", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1"), "Unknown flag -l");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "-n", "251"), "The specified flag \"--max-top-results\" has invalid value \"251\". It should be less than or equal to \"250\"");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "-n", "0"), "The specified flag \"--max-top-results\" has invalid value \"0\". It should be greater than or equal to \"1\"");
      assertStringContains(checkHandleFlagsErr("-t", "blaht", "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "13"), "The specified SDF, \"blaht\", does not exist");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-identity", "-1"), "Invalid minimum percentage identity -1 must be between 0 and 100");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-identity", "101"), "Invalid minimum percentage identity 101 must be between 0 and 100");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "101", "--min-bit-score", "5"), "Cannot set e-score and bit-score filters at the same time");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "-0.01", "--min-bit-score", "-1"), "Cannot set e-score and bit-score filters at the same time");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "-0.01"), "Invalid maximum e-score " + (-0.01) + " must be positive");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-bit-score", "-0.01"), "Invalid minimum bit score " + (-0.01) + " must be positive");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--min-bit-score", "0.5", "-n", "0"), "The specified flag \"--max-top-results\" has invalid value \"0\". It should be greater than or equal to \"1\"");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "12", "--max-e-score", "0.5", "-n", "0"), "The specified flag \"--max-top-results\" has invalid value \"0\". It should be greater than or equal to \"1\"");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "1", "-b", "1", "-w", "13"), "The specified flag \"--word\" has invalid value \"13\". It should be less than or equal to \"12\".");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "-1", "-b", "0", "-w", "12"), "The specified flag \"-a\" has invalid value \"-1\". It should be greater than or equal to \"0\".");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-a", "0", "-b", "-1", "-w", "12"), "The specified flag \"-b\" has invalid value \"-1\". It should be greater than or equal to \"0\".");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-o", output.getPath(), "-c", "0"), "The specified flag \"-c\" has invalid value \"0\". It should be greater than or equal to \"1\".");
      assertStringContains(checkHandleFlagsErr("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=blosum45", "--step", "-1"), "Unknown flag --step");

      checkMainInitOk("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=blosum45");
      checkMainInitOk("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", new File(dir, "output2").getPath(), "-w", "4", "-T", "1", "--matrix=blosum80");
      FileHelper.deleteAll(output);
      checkMainInitOk("-t", template.getPath(), "-i", left.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix=blosum45", "-f", "topn");

      final File right = new File(reads, "right");
      assertTrue(right.mkdir());
      ReaderTestUtils.getReaderDNA(READS_FASTA_ONE_INDEL, right, null).close();

      assertTrue(checkHandleFlagsErr("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", new File(dir, "output3").toString(), "-w", "4", "-T", "1", "--matrix=blosum45").contains("Paired end data not supported"));
    }
  }

}
