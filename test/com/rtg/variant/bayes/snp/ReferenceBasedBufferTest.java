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
package com.rtg.variant.bayes.snp;

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelFactory;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ReferenceBasedBufferTest extends TestCase {

  private static class MockModel extends Model<Description> {
    private final int mId;
    private final int mRef;

    MockModel(final int id, final int ref) {
      super(new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), true, ref), new StatisticsSnp(new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), true, ref).description()));
      mId = id;
      mRef = ref;
    }

    @Override
    public String toString() {
      return mId + ":" + mRef;
    }
  }

  private static class Fac implements ModelFactory<Description, Hypotheses<Description>> {
    private int mCount = 0;
    @Override
    public MockModel make(final int ref) {
      MockModel model;
      try {
        model = new MockModel(mCount++, ref);
      } catch (final InvalidParamsException e) {
        throw new RuntimeException();
      }
      return model;
    }
    @Override
    public Hypotheses<Description> defaultHypotheses(int ref) {
      throw new UnsupportedOperationException();
    }
  }

  private static final String INIT = ""
      + "Buffer length=1 base=0 current=0" + LS
      + ">>>>" + LS
      + "[0]     null" + LS
      + "" + LS
      ;
  private static final String EXP_3 = ""
      + "Buffer length=10 base=0 current=0" + LS
      + ">>>>" + LS
      + "[0]     null" + LS
      + "[1]     null" + LS
      + "[2]     null" + LS
      + "[3]     0:3" + LS
      + "[4]     null" + LS
      + "[5]     null" + LS
      + "[6]     null" + LS
      + "[7]     null" + LS
      + "[8]     null" + LS
      + "[9]     null" + LS
      + "" + LS
      ;
  private static final String EXP_SET = ""
      + "Buffer length=10 base=1 current=1" + LS
      + "[0]     null" + LS
      + ">>>>" + LS
      + "[1]     null" + LS
      + "[2]     null" + LS
      + "[3]     0:3" + LS
      + "[4]     null" + LS
      + "[5]     null" + LS
      + "[6]     null" + LS
      + "[7]     null" + LS
      + "[8]     null" + LS
      + "[9]     null" + LS
      + "" + LS
      ;
  private static final String EXP_SET1 = ""
      + "Buffer length=10 base=1 current=1" + LS
      + "[0]     null" + LS
      + ">>>>" + LS
      + "[1]     3:1" + LS
      + "[2]     null" + LS
      + "[3]     0:3" + LS
      + "[4]     null" + LS
      + "[5]     null" + LS
      + "[6]     null" + LS
      + "[7]     null" + LS
      + "[8]     4:-1" + LS
      + "[9]     2:3" + LS
      + "" + LS
      ;
  private static final String EXP_SET2 = ""
      + "Buffer length=24 base=1 current=0" + LS
      + ">>>>" + LS
      + "[0]     3:1" + LS
      + "[1]     null" + LS
      + "[2]     0:3" + LS
      + "[3]     null" + LS
      + "[4]     null" + LS
      + "[5]     null" + LS
      + "[6]     null" + LS
      + "[7]     4:-1" + LS
      + "[8]     2:3" + LS
      + "[9]     null" + LS
      + "[10]     null" + LS
      + "[11]     null" + LS
      + "[12]     null" + LS
      + "[13]     null" + LS
      + "[14]     null" + LS
      + "[15]     null" + LS
      + "[16]     null" + LS
      + "[17]     null" + LS
      + "[18]     null" + LS
      + "[19]     null" + LS
      + "[20]     null" + LS
      + "[21]     null" + LS
      + "[22]     null" + LS
      + "[23]     null" + LS
      + "" + LS
      ;
  private static final String EXP_STEP1 = ""
      + "Buffer length=24 base=24 current=23" + LS
      + "[0]     null" + LS
      + "[1]     null" + LS
      + "[2]     null" + LS
      + "[3]     null" + LS
      + "[4]     null" + LS
      + "[5]     null" + LS
      + "[6]     null" + LS
      + "[7]     null" + LS
      + "[8]     null" + LS
      + "[9]     null" + LS
      + "[10]     null" + LS
      + "[11]     null" + LS
      + "[12]     null" + LS
      + "[13]     null" + LS
      + "[14]     null" + LS
      + "[15]     null" + LS
      + "[16]     null" + LS
      + "[17]     null" + LS
      + "[18]     null" + LS
      + "[19]     null" + LS
      + "[20]     null" + LS
      + "[21]     null" + LS
      + "[22]     null" + LS
      + ">>>>" + LS
      + "[23]     null" + LS
      + "" + LS
      ;
  private static final String EXP_STEP2 = ""
      + "Buffer length=24 base=25 current=0" + LS
      + ">>>>" + LS
      + "[0]     null" + LS
      + "[1]     null" + LS
      + "[2]     null" + LS
      + "[3]     null" + LS
      + "[4]     null" + LS
      + "[5]     null" + LS
      + "[6]     null" + LS
      + "[7]     null" + LS
      + "[8]     null" + LS
      + "[9]     null" + LS
      + "[10]     null" + LS
      + "[11]     null" + LS
      + "[12]     null" + LS
      + "[13]     null" + LS
      + "[14]     null" + LS
      + "[15]     null" + LS
      + "[16]     null" + LS
      + "[17]     null" + LS
      + "[18]     null" + LS
      + "[19]     null" + LS
      + "[20]     null" + LS
      + "[21]     null" + LS
      + "[22]     null" + LS
      + "[23]     null" + LS
      + "" + LS
      ;

  public void test() {

    final byte[] template = {1, 2, 3, 4, 3, 2, 1, 2, 0, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1};
    final ReferenceBasedBuffer<ModelInterface<Description>> cb = new ReferenceBasedBuffer<>(1, new Fac(), template, 0);
    cb.globalIntegrity();
    assertEquals(INIT, cb.toString());
    assertEquals(0, cb.base());

    assertEquals("0:3", cb.get(3).toString());
    cb.globalIntegrity();

    assertEquals(EXP_3, cb.toString());
    assertEquals(0, cb.base());

    final ModelInterface<Description> step = cb.step();
    assertEquals(1, cb.base());

    assertEquals("1:0", step.toString());
    assertEquals(EXP_SET, cb.toString());
    assertEquals(1, cb.find(1));
    assertEquals(9, cb.find(9));
    assertEquals(0, cb.find(10));
    assertEquals(EXP_SET, cb.toString());
    assertEquals("2:3", cb.get(9).toString());
    assertEquals("3:1", cb.get(1).toString());
    assertNotNull(cb.get(8)); // N
    assertEquals(EXP_SET1, cb.toString());

    assertEquals(10, cb.find(11));
    cb.globalIntegrity();

    assertEquals(EXP_SET2, cb.toString());

    badFind(cb, -1, 1);
    badFind(cb, 0, 1);

    for (int i = 1; i < 24; i++) {
      final ModelInterface<Description> model = cb.step();
      if (template[i] >= 0) {
        assertNotNull("" + i, model);
      } else {
        assertNull("" + i, model);
      }
      cb.globalIntegrity();
    }
    assertEquals(EXP_STEP1, cb.toString());
    assertEquals(23, cb.find(24));
    assertEquals(0, cb.find(25));
    cb.step();
    cb.globalIntegrity();
    assertEquals(EXP_STEP2, cb.toString());
    badFind(cb, 1, 25);
  }

  private void badFind(final ReferenceBasedBuffer<?> mb, final int f, final int b) {
    try {
      mb.find(f);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Index less than base. index=" + Integer.toString(f) + " base=" + Integer.toString(b), e.getMessage());
    }
  }
}
