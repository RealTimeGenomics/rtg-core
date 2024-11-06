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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.Description;

import junit.framework.TestCase;

/**
 */
public class DescriptionSnpTest extends TestCase {

  private void check(final String exp, final int index, Description description) {
    assertTrue(description.valid(index));
    assertEquals(exp, description.name(index));
    final StringBuilder sb = new StringBuilder();
    description.writeName(sb, index);
    assertEquals(exp, sb.toString());
  }

  public void test() {
    assertEquals(4, DescriptionSnp.SINGLETON.size());
    assertEquals(1, DescriptionSnp.SINGLETON.minLength());
    assertEquals(1, DescriptionSnp.SINGLETON.maxLength());
    check("A", 0, DescriptionSnp.SINGLETON);
    check("C", 1, DescriptionSnp.SINGLETON);
    check("G", 2, DescriptionSnp.SINGLETON);
    check("T", 3, DescriptionSnp.SINGLETON);

    assertFalse(DescriptionSnp.SINGLETON.valid(-1));
    assertFalse(DescriptionSnp.SINGLETON.valid(4));
  }
}
