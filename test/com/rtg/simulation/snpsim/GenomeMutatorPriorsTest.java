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
package com.rtg.simulation.snpsim;

import java.io.File;

import com.rtg.reference.Ploidy;
import com.rtg.simulation.snpsim.Mutation.DifferentMode;
import com.rtg.simulation.snpsim.Mutation.MutationType;
import com.rtg.util.PortableRandom;
import com.rtg.util.TestUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;

import junit.framework.TestCase;

/**
 */
public class GenomeMutatorPriorsTest extends TestCase {

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(GenomeMutatorPriorsTest.class);
  }
  private GenomePriorParams mPriors;
  private File mDir;

  @Override
  public void setUp() throws Exception {
    mPriors = GenomePriorParams.builder().genomePriors("testhumanprior").create();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    mPriors = null;
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  public void testSimpleRates() {
    final GenomeMutatorPriors mg = new GenomeMutatorPriors(0.5, 0.0, 0.5, 0.1);
    final String str = mg.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, "SnpThres " + 0.5, "MnpThres " + 0.5,
              "InsertThres " + 0.75, "DeleteThres " + 1.0, "InsDelThres " + 1.0,
              "SnpEteroThres " + 0.33, "MnpEteroThres " + 0.33, "InsertEteroThres " + 0.33,
              "DeleteEteroThres " + 0.33);
  }

  public void testOnlyIndels() {
    final GenomeMutatorPriors mg = new GenomeMutatorPriors(0.0, 0.0, 1.0, 0.1);
    final String str = mg.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, "SnpThres " + 0.0, "MnpThres " + 0.0,
              "InsertThres " + 0.5, "DeleteThres " + 1.0, "InsDelThres " + 1.0,
              "SnpEteroThres " + 0.33, "MnpEteroThres " + 0.33, "InsertEteroThres " + 0.33,
              "DeleteEteroThres " + 0.33);
    MutationType t = mg.chooseType(0.2);
    assertTrue(t == MutationType.INSERT);
    t = mg.chooseType(0.6);
    assertTrue(t == MutationType.DELETE);
  }

  private void checkAccumDist(final GenomeMutatorPriors mg, final MutationType type) {
    double r = 0.0;
    double last = 0;
    for (int i = 0; i < 22; i++) {
      if (r > 1.0) {
        r = 1.0;
      }
      final int len = mg.chooseLength(r, type, true);
      //System.err.println("rand MNP " + r + " len " + len);
      assertTrue(len > 0);
      assertTrue(len <= mg.maxLength(type, true));
      assertTrue(len >= last);
      last = len;
      r += 0.05;
    }
  }

  public void testPriorRates() {
    final GenomeMutatorPriors mg = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    final String str = mg.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, "SnpEteroThres " + 0.5);
    final DifferentMode mode = mg.chooseDifferentMode(0.2, MutationType.SNP);
    assertTrue(mode == DifferentMode.HOMOZYGOUS);
    final DifferentMode mode2 = mg.chooseDifferentMode(0.6, MutationType.SNP);
    assertFalse(mode2 == DifferentMode.HOMOZYGOUS);
    checkAccumDist(mg, MutationType.DELETE);
    checkAccumDist(mg, MutationType.MNP);
    final double[][] snpProbabilities = mg.snpProbabilities();
    assertEquals(4, snpProbabilities.length);
    for (double[] snpProbability : snpProbabilities) {
      assertEquals(9, snpProbability.length);
      for (final double p : snpProbability) {
        assertTrue(p > 0);
        assertTrue(p < 0.1);
      }
    }
    byte[] snps = mg.chooseAltSnp(0.0, (byte) 0);
    assertEquals(2, snps.length);
    assertEquals(0, snps[0]);
    assertEquals(0, snps[1]);
    for (byte i = 1; i <= 4; i++) {
      try {
        mg.chooseAltSnp(1.1, i);
        fail();
      } catch (final IllegalArgumentException e) {
        assertEquals("Invalid snp distribution", e.getMessage());
      }
      snps = mg.chooseAltSnp(0.99999999, i);
      assertEquals(2, snps.length);
      assertEquals(3, snps[0]);
      assertEquals(4, snps[1]);
    }


  }


  public void testDefaultIndelLengths() {
    final GenomeMutatorPriors priors = new GenomeMutatorPriors(0.2, 0.2, 0.2, 0.6);
    final PortableRandom r = new PortableRandom(123);
    final int[] counts = new int[GenomeMutatorPriors.DEFAULT_INDEL_DIST.length];
    double total = 0;
    for (int i = 0; i < counts.length; i++) {
      counts[i] = 0;
      total += GenomeMutatorPriors.DEFAULT_INDEL_DIST[i];
    }

    assertTrue(Math.abs(total - 1) < 0.000001);


    final int iterations = 100000;
    for (int i = 0; i < iterations; i++) {
      counts[priors.chooseLength(r, MutationType.INSERT, false) - 1]++;
    }

    assertEquals(0, counts[0]);

    assertTrue(counts[1] < iterations * 0.48 * 1.1);
    assertTrue(counts[1] > iterations * 0.48 * 0.9);

    assertTrue(counts[2] < iterations * 0.17 * 1.1);
    assertTrue(counts[2] > iterations * 0.17 * 0.9);

    assertTrue(counts[3] < iterations * 0.05 * 1.1);
    assertTrue(counts[3] > iterations * 0.05 * 0.9);

    assertTrue(counts[4] < iterations * 0.09 * 1.1);
    assertTrue(counts[4] > iterations * 0.09 * 0.9);

    assertTrue(counts[5] < iterations * 0.02 * 1.1);
    assertTrue(counts[5] > iterations * 0.02 * 0.9);

    assertTrue(counts[counts.length - 1] < iterations * 0.015 * 1.1);
    assertTrue(counts[counts.length - 1] > iterations * 0.015 * 0.9);
  }

  public void testToString() {

    GenomeMutatorPriors priors = new GenomeMutatorPriors(0.2, 0.3, 0.5, 1.0);

    String result = priors.toString();
    TestUtils.containsAll(result,
              "SnpThres " + Double.toString(0.2),
              "MnpThres " + Double.toString(0.5),
              "InsertThres " + Double.toString(0.75),
              "InsDelThres " + Double.toString(1),
              "SnpEteroThres " + Double.toString(0.33),
              "MnpEteroThres " + Double.toString(0.33),
              "InsertEteroThres " + Double.toString(0.33),
              "DeleteEteroThres " + Double.toString(0.33),
              "VariantDifferentFactor " + Double.toString(0.33333333333),
              "InsDelEteroThres n/a",
              "MnpOmoDist",
              "MnpEteroDist",
              "IndelOmoDist",
              "IndelEteroDist"
    );

    priors = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    result = priors.toString();
    //System.err.println(result);
    TestUtils.containsAll(result,
              "SnpThres 0.9018" ,
              "MnpThres 0.9368" ,
              "InsertThres " ,
              "InsDelThres " ,
              "SnpEteroThres " ,
              "MnpEteroThres " ,
              "InsertEteroThres " ,
              "DeleteEteroThres " ,
              "VariantDifferentFactor " + Double.toString(0.33333333333),
              "InsDelEteroThres n/a",
              "MnpOmoDist",
              "MnpEteroDist",
              "IndelOmoDist",
              "IndelEteroDist",
              "MnpThresholds-ends: [" + Double.toString(0.5) + ", " + Double.toString(0.75) + ", " + Double.toString(1.0) + ", " + Double.toString(1.0) + "]",
              "MnpThresholds-ends-snponly: [" + Double.toString(1.0) + ", " + Double.toString(1.0) + ", " + Double.toString(1.0) + ", " + Double.toString(1.0) + "]",
              "MnpThresholds-mid1: [" + Double.toString(0.05) + ", " + Double.toString(0.1) + ", " + Double.toString(0.15) + ", " + Double.toString(1.0) + "]",
              "MnpThresholds-mid2: [" + Double.toString(0.5) + ", " + Double.toString(0.7) + ", " + Double.toString(0.9) + ", " + Double.toString(1.0) + "]"
    );


    //assertTrue(result.contains("InsertThres 0.5"));
  }

  private void validateTypeRates(final PortableRandom r, final GenomeMutatorPriors priors, final double[] rates) {

    final int iterations = 100000;
    int mnp = 0;
    int snp = 0;
    int indel = 0;
    for (int i = 0; i < iterations; i++) {
      final MutationType type = priors.chooseType(r);
      if (type == MutationType.SNP) {
        snp++;
      } else if (type == MutationType.MNP) {
        mnp++;
      } else {
        indel++;
      }
    }
    assertTrue(Math.abs(snp - iterations * rates[0]) <= snp * 0.1);
    assertTrue(Math.abs(mnp - iterations * rates[1]) <= mnp * 0.1);
    assertTrue(Math.abs(indel - iterations * rates[2]) <= indel * 0.1);

  }

  private void validateDiffModeRates(final PortableRandom r, final MutationType type, final GenomeMutatorPriors priors, final double[] rates) {
    final int iterations = 100000;
    int omo = 0;
    int one = 0;
    int diff = 0;
    for (int i = 0; i < iterations; i++) {
       final DifferentMode mode = priors.chooseDifferentMode(r, type);

       switch (mode) {
         case HOMOZYGOUS:
           omo++;
           break;
         case ONE_ONLY:
           one++;
           break;
         case DIFFERENT:
           diff++;
           break;
         default:
           throw new RuntimeException();
       }
    }
    assertTrue(Math.abs(omo - iterations * rates[0]) <= omo * 0.1);
    assertTrue(Math.abs(one - iterations * rates[1]) <= one * 0.1);
    assertTrue(Math.abs(diff - iterations * rates[2]) <= diff * 0.1);

  }

  public void testDefaultRates() {
    final PortableRandom r = new PortableRandom(123);
    final double[] rates = {0.2, 0.3, 0.5};
    final double third = 1.0 / 3;
    final double[] thirdRates = {third, third, third};

    GenomeMutatorPriors priors = new GenomeMutatorPriors(0.2, 0.3, null, 1.0);
    validateTypeRates(r, priors, rates);

    priors = new GenomeMutatorPriors(0.2, null, 0.5, 1.0);
    validateTypeRates(r, priors, rates);

    priors = new GenomeMutatorPriors(null, 0.3, 0.5, 1.0);
    validateTypeRates(r, priors, rates);

    priors = new GenomeMutatorPriors(0.0, 0.0, 0.0, 1.0);
    validateTypeRates(r, priors, thirdRates);

    priors = new GenomeMutatorPriors(0.0, 0.1, 0.4, 1.0);
    validateTypeRates(r, priors, new double[] {0.0, 0.2, 0.8});
  }

  public void testPriorParamsConstructor() {
    final PortableRandom r = new PortableRandom(123);
    final GenomeMutatorPriors priors = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    final double snpRate = 0.901826;
    final double mnpRate = 0.035046;
    final double[] rates = {snpRate, mnpRate, 1 - snpRate - mnpRate};
    validateTypeRates(r, priors, rates);

    final double[] snpMode = {0.5, 0.5 * 0.66, 0.5 * 0.33};
    validateDiffModeRates(r, MutationType.SNP, priors, snpMode);

    final double[] mnpMode = {0.38, 0.622 * 0.66, 0.622 * 0.33};
    validateDiffModeRates(r, MutationType.MNP, priors, mnpMode);


    final double[] insertMode = {0.6996, 0.24031, 0.06005};
    validateDiffModeRates(r, MutationType.INSERT, priors, insertMode);

    final double[] delMode = {0.6996, 0.24031, 0.06005};
    validateDiffModeRates(r, MutationType.DELETE, priors, delMode);

    final double[] insdelMode = {0, 0, 1};
    validateDiffModeRates(r, MutationType.INSDEL, priors, insdelMode);

    assertTrue(Math.abs(priors.rate() - 0.00066531) < 0.00000001);

  }

  public void testMaxLength() {
    final GenomeMutatorPriors priors = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    assertEquals(11, priors.maxLength(MutationType.INSERT, true));
    assertEquals(11, priors.maxLength(MutationType.INSERT, false));
    assertEquals(11, priors.maxLength(MutationType.DELETE, true));
    assertEquals(11, priors.maxLength(MutationType.DELETE, false));
    assertEquals(7, priors.maxLength(MutationType.MNP, true));
    assertEquals(7, priors.maxLength(MutationType.MNP, false));
    try {
      assertEquals(1, priors.maxLength(MutationType.SNP, false));
      fail();
    } catch (final Exception e) {
      assertEquals("Unpossible", e.getMessage());
    }
    try {
      assertEquals(1, priors.maxLength(MutationType.SNP, true));
      fail();
    } catch (final Exception e) {
      //expected
      assertEquals("Unpossible", e.getMessage());
    }

  }
}
