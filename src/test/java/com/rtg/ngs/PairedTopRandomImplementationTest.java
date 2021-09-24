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

package com.rtg.ngs;

import com.rtg.ngs.PairedTopRandomImplementation.HitRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class PairedTopRandomImplementationTest extends TestCase {

  public PairedTopRandomImplementation getImpl() {
    return new PairedTopRandomImplementation(3, 7);
  }

  public void testMost() {
    final PairedTopRandomImplementation petri = getImpl();
    assertNull(petri.getRecords());

    petri.initialize();

    assertEquals(3, petri.getRecords().length);

    petri.update(1, 0, 1, true, 1, false, 5, 1);

    HitRecord hr = petri.getRecords()[0];
    assertNull(hr);

    hr = petri.getRecords()[1];
    assertNotNull(hr);
    assertEquals(0, hr.mTemplateId);
    assertEquals(1, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(1, hr.mZeroBasedMateTemplateStart);
    assertFalse(hr.mMateReverse);
    assertEquals(5, hr.mComboScore);

    hr = petri.getRecords()[2];
    assertNull(hr);

    petri.update(2, 0, 4, true, 1, true, 6, 1);  //worse score should not overwrite
    hr = petri.getRecords()[1];
    assertNotNull(hr);
    assertEquals(0, hr.mTemplateId);
    assertEquals(1, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

  }

  public void testRandomness() {
    final PairedTopRandomImplementation petri = getImpl();
    petri.initialize();
    petri.update(0, 0, 1, true, 1, true, 5, 1);

    HitRecord hr = petri.getRecords()[0];
    assertNotNull(hr);
    assertEquals(0, hr.mTemplateId);
    assertEquals(1, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

    petri.update(0, 1, 2, true, 1, true, 5, 1);  //n=1 definitely overwrites
    hr = petri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

    petri.update(0, 2, 3, true, 1, true, 5, 2);  //n=2 50/50 overwrites, with this seed this call does not
    hr = petri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

    petri.update(0, 2, 3, true, 1, true, 5, 2);  //try again, same
    hr = petri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

    petri.update(0, 2, 3, true, 1, true, 5, 2);  //try again, same
    hr = petri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

    petri.update(0, 2, 3, true, 1, true, 5, 2);  //try again, successfully overwritten
    hr = petri.getRecords()[0];
    assertEquals(2, hr.mTemplateId);
    assertEquals(3, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mComboScore);

  }
}
