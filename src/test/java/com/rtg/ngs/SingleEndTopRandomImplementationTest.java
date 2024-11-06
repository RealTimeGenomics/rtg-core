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

package com.rtg.ngs;

import com.rtg.ngs.SingleEndTopRandomImplementation.HitRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SingleEndTopRandomImplementationTest extends TestCase {

  public SingleEndTopRandomImplementation getImpl() {
    return new SingleEndTopRandomImplementation(3, 7);
  }

  public void testMost() {
    final SingleEndTopRandomImplementation setri = getImpl();
    assertNull(setri.getRecords());

    setri.initialize();

    assertEquals(3, setri.getRecords().length);

    setri.update(1, 0, 1, true, 5, 1);

    HitRecord hr = setri.getRecords()[0];
    assertNull(hr);

    hr = setri.getRecords()[1];
    assertNotNull(hr);
    assertEquals(0, hr.mTemplateId);
    assertEquals(1, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

    hr = setri.getRecords()[2];
    assertNull(hr);

    setri.update(2, 0, 4, true, 6, 1);  //worse score should not overwrite
    hr = setri.getRecords()[1];
    assertNotNull(hr);
    assertEquals(0, hr.mTemplateId);
    assertEquals(1, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

  }

  public void testRandomness() {
    final SingleEndTopRandomImplementation setri = getImpl();
    setri.initialize();
    setri.update(0, 0, 1, true, 5, 1);

    HitRecord hr = setri.getRecords()[0];
    assertNotNull(hr);
    assertEquals(0, hr.mTemplateId);
    assertEquals(1, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

    setri.update(0, 1, 2, true, 5, 1);  //n=1 definitely overwrites
    hr = setri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

    setri.update(0, 2, 3, true, 5, 2);  //n=2 50/50 overwrites, with this seed this call does not
    hr = setri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

    setri.update(0, 2, 3, true, 5, 2);  //try again, same
    hr = setri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

    setri.update(0, 2, 3, true, 5, 2);  //try again, same
    hr = setri.getRecords()[0];
    assertEquals(1, hr.mTemplateId);
    assertEquals(2, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);

    setri.update(0, 2, 3, true, 5, 2);  //try again, successfully overwritten
    hr = setri.getRecords()[0];
    assertEquals(2, hr.mTemplateId);
    assertEquals(3, hr.mZeroBasedTemplateStart);
    assertTrue(hr.mReverse);
    assertEquals(5, hr.mAlignScore);
  }
}
