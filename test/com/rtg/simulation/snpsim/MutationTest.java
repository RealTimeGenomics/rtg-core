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

import com.rtg.reference.Ploidy;
import com.rtg.util.PortableRandom;
import com.rtg.util.TestUtils;
import com.rtg.variant.GenomePriorParams;

import junit.framework.TestCase;

/**
 */
public class MutationTest extends TestCase {

  private GenomePriorParams mPriors;

  @Override
  public void setUp() throws Exception {
    mPriors = GenomePriorParams.builder().genomePriors("testhumanprior").create();
  }

  @Override
  public void tearDown() {
    mPriors = null;
  }

  static boolean integrity(final Mutation m) {
    if (!m.mHeterozygous) {
      // is homozygous
      if ((m.mDiffMode != Mutation.DifferentMode.HOMOZYGOUS) || (m.mGenDiffMode != GenDiffMode.BOTH_SAME)) {
        return false;
      }
      if (m.mGenDiffMode == GenDiffMode.BOTH_SAME) {
        if (m.mLength != m.mLengthTwin) {
          return false;
        }
      }
      return true;
    } else {
      // is heterozygous
      if ((m.mDiffMode == Mutation.DifferentMode.HOMOZYGOUS) || (m.mGenDiffMode == GenDiffMode.BOTH_SAME)) {
        return false;
      }

      if (m.mType == Mutation.MutationType.DELETE && m.mLength == m.mLengthTwin) {
        return false;
      }
      if (m.mGenDiffMode == GenDiffMode.DIFFERENT) {
        if (m.mLength == 0 || m.mLengthTwin == 0) {
          return false;
        }
      }
      // SNPs are both sides length 1 even if only one side is changed
      // MNPs the unchanged side has the length of the effected bases for the
      // other side
      if ((m.mType == Mutation.MutationType.INSERT) || (m.mType == Mutation.MutationType.DELETE) || (m.mType == Mutation.MutationType.INSDEL)) {
        if (m.mGenDiffMode == GenDiffMode.FIRST_ONLY) {
          if (m.mLength == 0 || m.mLengthTwin != 0) {
            return false;
          }
        }
        if (m.mGenDiffMode == GenDiffMode.TWIN_ONLY) {
          if (m.mLength != 0 || m.mLengthTwin == 0) {
            return false;
          }
        }
      }
    }
    return true;
  }

  public void test() {
    final Mutation m1 = new Mutation(1, Mutation.MutationType.SNP, true,
        Mutation.DifferentMode.ONE_ONLY, GenDiffMode.FIRST_ONLY, 2, 0,
      null, null);
    assertTrue(integrity(m1));
    final GenomeMutatorPriors mg = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    final Mutation m2 = new Mutation(new PortableRandom(1), 1, mg, 200, (byte) 0, null);
    assertTrue(integrity(m2));
    final String str = m2.makeString();
    assertTrue(str, str.contains("   pos:"));
  }

  public void testMakeString() {
    final Mutation m1 = new Mutation(20, Mutation.MutationType.SNP, true, Mutation.DifferentMode.DIFFERENT, GenDiffMode.DIFFERENT, 1, 2, null, null);
    assertEquals("   pos:20 e SNP l1: 1 l2: 2 DIFFERENT", m1.makeString());
    assertTrue(integrity(m1));
  }

  public void testEnum() {
    TestUtils.testEnum(Mutation.DifferentMode.class, "[HOMOZYGOUS, ONE_ONLY, DIFFERENT]");
    TestUtils.testPseudoEnum(Mutation.MutationType.class, "[SNP, MNP, INSERT, DELETE, INSDEL, PREPARED]");
  }

  private int sum(final byte[] n) {
    int s = 171; // just so null's contribute
    if (n != null) {
      for (final byte b : n) {
        s += b;
      }
    }
    return s;
  }

  public void testRandomGeneration() {
    final GenomeMutatorPriors mg = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    // Use a seed to avoid intermittent failure if thresholds aren't quite wide enough.
    final PortableRandom r = new PortableRandom(500);
    int snpCount = 0;

    int mnpCount = 0;
    final int iterations = 100000;

    int heteroCount = 0;

    int firstDel = 0;
    int notFirstDel = 0;

    int twin = 0;
    int first = 0;

    int sum = 0;

    for (int i = 0; i < iterations; i++) {
      final Mutation m = new Mutation(r, 1000, mg, 2000, (byte) (i % 4 + 1), null);
      assertTrue(integrity(m));
      sum += sum(m.mBases);
      if (m.mDiffMode == Mutation.DifferentMode.ONE_ONLY) {
        assertTrue(m.mGenDiffMode != GenDiffMode.DIFFERENT);
      }
      if (m.mType == Mutation.MutationType.SNP) {
        snpCount++;
        assertEquals(1, m.mLength);
        assertEquals(1, m.mLengthTwin);
      } else if (m.mType == Mutation.MutationType.DELETE) {
        if (m.mDiffMode != Mutation.DifferentMode.HOMOZYGOUS) {
          assertTrue(m.mLength != m.mLengthTwin);
        }
      } else if (m.mType == Mutation.MutationType.MNP) {
        mnpCount++;
      }

      if (m.mType != Mutation.MutationType.SNP) {
        if (m.isFirstDelete()) {
          firstDel++;
        } else {
          notFirstDel++;
        }
      } else {
        // has no meaning so should be true
        assertTrue(m.isFirstDelete());
      }


      if (m.mHeterozygous) {
        heteroCount++;
        assertTrue(m.makeString().startsWith("   pos:1000 e"));
        if (m.mGenDiffMode == GenDiffMode.FIRST_ONLY) {
          first++;
          if (m.mType != Mutation.MutationType.SNP) {
            assertEquals(0, m.mLengthTwin);
          }
        } else if (m.mGenDiffMode == GenDiffMode.TWIN_ONLY) {
          twin++;
          if (m.mType != Mutation.MutationType.SNP) {
            assertEquals(0, m.mLength);
          }
        }
      } else {
        assertTrue(m.mLength == m.mLengthTwin);
        assertTrue(m.makeString().startsWith("   pos:1000 o"));
        assertEquals(GenDiffMode.BOTH_SAME, m.mGenDiffMode);
      }
    }

    assertEquals(17550783, sum); // REGRESSION

    //Calculated manually from the testhumanpriors file
    // Indels will make up the remainder but are split into 3 types.

    // testhuman priors gives 0.90183 snps 0.0003*2/(0.0003 * 2 + 0.00008 / 5.8 + 0.00002 / 2.1 + 0.000042)
    assertTrue(snpCount > 0.90183 * iterations * 0.9);
    assertTrue(snpCount < 0.90183 * iterations * 1.1);

    // 0.03504 = (0.00008 / 5.8 + 0.00002 / 2.1)/(0.0003 * 2 + 0.00008 / 5.8 + 0.00002 / 2.1 + 0.000042)
    assertTrue(mnpCount > 0.03504 * iterations * 0.9);
    assertTrue(mnpCount < 0.03504 * iterations * 1.1);

    // 0.59833 = (0.00008 / 5.8 + 0.000807 / 2.1)/(0.0003 * 2 + 0.00008 / 5.8 + 0.00002 / 2.1 + 0.000042)
    assertTrue(heteroCount > 0.59833 * iterations * 0.9);
    assertTrue(heteroCount < 0.59833 * iterations * 1.1);

    // these should be roughly equal
    assertTrue(notFirstDel > firstDel - (firstDel / 20));
    assertTrue(notFirstDel < firstDel + (firstDel / 20));

    //twin only should be ~= first
    assertTrue(twin > first - (first / 20));
    assertTrue(twin < first + (first / 20));

  }

  public void testRandomEndClipping() {
    final GenomeMutatorPriors mg = new GenomeMutatorPriors(mPriors, Ploidy.DIPLOID);
    final PortableRandom r = new PortableRandom();
    for (int i = 0; i < 50000; i++) {
      final Mutation m = new Mutation(r, 1000, mg, 1002, (byte) 0, null);
      assertTrue(integrity(m));
      if (m.mType != Mutation.MutationType.INSERT) {
        assertTrue(m.mLength <= 2);
        assertTrue(m.mLengthTwin <= 2);
      }
    }
    for (int i = 0; i < 50000; i++) {
      final Mutation m = new Mutation(r, 1001, mg, 1002, (byte) 0, null);
      assertTrue(integrity(m));
      if (m.mType != Mutation.MutationType.INSERT) {
        assertTrue(m.mLength <= 1);
        assertTrue(m.mLengthTwin <= 1);
      }
      if (m.mType == Mutation.MutationType.DELETE) {
        assertTrue(m.mDiffMode == Mutation.DifferentMode.HOMOZYGOUS || m.mDiffMode == Mutation.DifferentMode.ONE_ONLY);
      }
    }
  }

}
