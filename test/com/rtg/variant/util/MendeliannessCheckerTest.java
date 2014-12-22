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
package com.rtg.variant.util;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class MendeliannessCheckerTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new MendeliannessChecker();
  }

  public void testHelp() {
    checkHelp("rtg mendelian",
        "output only consistent calls",
        "output only non-Mendelian calls",
        "all records, regardless",
        "VCF file"
        );
  }

  static boolean isBadHaploidChild(final String fatherCall, final String motherCall, final String childCall) {
    return MendeliannessAnnotator.isBadHaploidChild(new Genotype(fatherCall), new Genotype(motherCall), new Genotype(childCall));
  }

  static boolean isBadDiploidChild(final String fatherCall, final String motherCall, final String childCall) {
    return MendeliannessAnnotator.isBadDiploidChild(new Genotype(fatherCall), new Genotype(motherCall), new Genotype(childCall));
  }

  public void testBadHaploidChild() {
    final String ref = "0";
    final String var = "1";
    final String missing = ".";
    final String missing2 = "./.";
    final String homref = "0/0";
    final String homvar = "1/1";
    final String hetvar = "0/1";

    assertFalse(isBadHaploidChild(ref, hetvar, var));  // E.g. Chr X
    assertFalse(isBadHaploidChild(homvar, ref, var));  // Just for symmetry
    assertFalse(isBadHaploidChild(hetvar, var, ref));  // Just for symmetry
    assertFalse(isBadHaploidChild(var, missing, var)); // E.g. Chr Y
    assertFalse(isBadHaploidChild(missing, var, var)); // Just for symmetry

    assertTrue(isBadHaploidChild(ref, homvar, ref));  // E.g. Chr X
    assertTrue(isBadHaploidChild(var, var, ref)); // E.g. Chr Y
    assertTrue(isBadHaploidChild(var, missing, ref));  // E.g. Chr Y

    // It would be acceptable to change the outcomes of these ones:
    assertTrue(isBadHaploidChild(homvar, homref, var));  // Edge case, both parents are diploid, but we prioritize mother
    assertTrue(isBadHaploidChild(homref, homref, var)); // Edge case, both parents are diploid
    assertTrue(isBadHaploidChild(homref, missing2, var)); // Edge case, both parents are diploid
  }

  public void testBadDiploidChild() {
    final String homref = "0/0";
    final String homvar = "1/1";
    final String hetvar = "0/1";
    assertFalse(isBadDiploidChild(hetvar, homref, hetvar));
    assertFalse(isBadDiploidChild(hetvar, hetvar, hetvar));
    assertFalse(isBadDiploidChild(hetvar, hetvar, homvar));

    assertTrue(isBadDiploidChild(homref, homref, hetvar));
    assertTrue(isBadDiploidChild(homref, homvar, homvar));
    assertTrue(isBadDiploidChild(homref, homvar, homref));
  }

  public void testOptions() throws IOException {
    try (TestDirectory dir = new TestDirectory("mendelianness")) {
      final File sdf = ReaderTestUtils.getDNADir(">chr21\nacgt", dir);
      final File file1 = FileHelper.resourceToFile("com/rtg/variant/bayes/multisample/family/resources/merge.vcf", new File(dir, "merge.vcf"));
      final File inconsistent = new File(dir, "failed.vcf");
      final File consistent = new File(dir, "nonfailed.vcf");
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        final MendeliannessChecker chk = (MendeliannessChecker) getCli();

        int ret = chk.mainInit(new String[] {"-t", sdf.getPath(), "-i", file1.getPath(), "-Z", "--all-records", "--output-inconsistent", inconsistent.getPath(), "--output-consistent", consistent.getPath()}, bos, new PrintStream(bos));
        assertEquals(0, ret);
      } finally {
        bos.close();
      }
      final String s = bos.toString().replaceAll("Checking: [^\n]*\n", "Checking: \n");
      //System.err.println(s);
      mNano.check("mendelian1", s);
      final String s2 = IOUtils.readAll(inconsistent);
      //System.err.println(s2);
      mNano.check("mendelian2", s2);
      final String s2b = IOUtils.readAll(consistent);
      //System.err.println(s2);
      mNano.check("mendelian2b", s2b);

      final ByteArrayOutputStream bos3 = new ByteArrayOutputStream();
      try {
        final MendeliannessChecker chk = (MendeliannessChecker) getCli();
        int ret = chk.mainInit(new String[] {"-t", sdf.getPath(), "-i", file1.getPath()}, bos3, new PrintStream(bos));
        assertEquals(0, ret);
      } finally {
        bos3.close();
      }
      final String s3 = bos3.toString().replaceAll("Checking: [^\n]*\n", "Checking: \n");
      //System.err.println(s3);
      mNano.check("mendelian3", s3);
    }

  }

}
