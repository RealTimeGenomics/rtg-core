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
package com.rtg.vcf;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.Resources;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.util.BlockCompressedInputStream;

/**
 * Test the corresponding class.
 */
public class VcfDecomposerCliTest extends AbstractCliTest {

  private static final String RESOURCES = "com/rtg/vcf/resources/";

  @Override
  protected AbstractCli getCli() {
    return new VcfDecomposerCli();
  }

  public void testFlags() {
    checkHelp("rtg vcfdecompose", "Decomposes complex variants within a VCF file."
      , "-i,", "--input=FILE", "VCF file containing variants"
      , "-o,", "--output=FILE", "output VCF file"
      , "-t", "--template=SDF"
      , "--no-gzip"
      , "--no-header"
      , "--no-index"
    );
  }

  // Run a test where vcfdecompose is expected to complete normally
  private void runResourceTest(String inResourceLoc, String expResourceLoc, String... extraArgs) throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File in = new File(dir, new File(Resources.getResource(inResourceLoc).getFile()).getName());

      try (InputStream stream = Resources.getResourceAsStream(inResourceLoc)) {
        assertNotNull("Cant find:" + inResourceLoc, stream);

        final String snps = FileUtils.streamToString(stream);
        FileUtils.stringToFile(snps, in);
      }
      final File out = new File(dir, "out.vcf.gz");
      final String output = checkMainInitOk(Utils.append(extraArgs, "-i", in.getPath(), "-o", out.getPath()));
      mNano.check(expResourceLoc + ".txt", output, true);

      assertEquals(BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK, BlockCompressedInputStream.checkTermination(out));

      final String o = StringUtils.grep(FileHelper.gzFileToString(out), "^[^#]").replaceAll("[\r\n]+", "\n");
      mNano.check(expResourceLoc, o, true);
    }
  }

  private static final String REF = ">1\n"
    + "TGGGGAAGCAAGGCGGAGTTGGGCAGCTCGTGTTCAATGGGTAGAGTTTCAGGCTGGGGTGATGGAAGGGTGCTGGAAATGAGTGGTAGTGATGGCGGCACAACAGTGTGAATCT"
    + "ACTTAATCCCACTGAACTGTATGCTGAAAAATGGTTTAGACGGTGAATTTTAGGTTATGTATGTTTTACCACAATTTTTAAAAAGCTAGTGAAAAGCTGGTAAAAAGAAAGAAAA"
    + "GAGGCTTTTTTAAAAAGTTAAATATATAAAAAGAGCATCATCAGTCCAAAGTCCAGCAGTTGTCCCTCCTGGAATCCGTTGGCTTGCCTCCGGCATTTTTGGCCCTTGCCTTTTA"
    + "GGGTTGCCAGATTAAAAGACAGGATGCCCAGCTAGTTTGAATTTTAGATAAACAACGAATAATTTCGTAGCATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAA"
    + "AACATTATTGGTTGTTTATCTGAGATTCAGAATTAAGCATTTTATATTTTATTTGCTGCCTCTGGCCACCCTACTCTCTTCCTAACACTCTCTCCCTCTCCCAGTTTTGTCCGCC"
    + "TTCCCTGCCTCCTCTTCTGGGGGAGTTAGATCGAGTTGTAACAAGAACATGCCACTGTCTCGCTGGCTGCAGCGTGTGGTCCCCTTACCAGAGGTAAAGAAGAGATGGATCTCCA"
    + "CTCATGTTGTAGACAGAATGTTTATGTCCTCTCCAAATGCTTATGTTGAAACCCTAACCCCTAATGTGATGGTATGTGGAGATGGGCCTTTGGTAGGTAATTACGGTTAGATGAG"
    + "GTCATGGGGTGGGGCCCTCATTATAGATCTGGTAAGAAAAGAGAGCATTGTCTCTGTGTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCATTTCTCTCTATCTCATTTCTCTCTCT"
    + "CTCGCTATCTCATTTTTCTCTCTCTCTCTTTCTCTCCTCTGTCTTTTCCCACCAAGTGAGGATGCGAAGAGAAGGTGGCT";

  public void testDecompose() throws IOException {
    final File sdf = ReaderTestUtils.getDNADir(REF);
    try {
      runResourceTest(RESOURCES + "test.vcf", "expected.vcf", "-t", sdf.getPath());
    } finally {
      assertTrue(FileUtils.deleteFiles(sdf));
    }
  }
}
