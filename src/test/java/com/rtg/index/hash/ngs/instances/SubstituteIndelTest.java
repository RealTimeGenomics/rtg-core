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
package com.rtg.index.hash.ngs.instances;



import java.io.IOException;
import java.util.Collection;

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.util.StringUtils;


/**
 */
public final class SubstituteIndelTest extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return null;
  }

  /**
   * Test the substitution testing code.
   */
  private void checkSubstitute(final int n, final int expected, final String expStr, final int maxIndel) {
    final String testString = "acg";
    final SubstituteIndel sub = new SubstituteIndel(testString, maxIndel);
    final Collection<String> coll = sub.substitute(n);
    final StringBuilder sb = new StringBuilder();
    for (final String templ : coll) {
      sb.append(templ).append(StringUtils.LS);
    }
    if (expStr != null) {
      assertEquals(expStr, sb.toString());
    }
    final int actualCount = coll.size();
    assertEquals(expected, actualCount);
  }

  private static final String EXPECTED_INDEL = ""
    + "cgt" + StringUtils.LS
    + "cacg" + StringUtils.LS
    + "agt" + StringUtils.LS
    + "accg" + StringUtils.LS
    + "aag" + StringUtils.LS
    + "act" + StringUtils.LS
    + "accg" + StringUtils.LS
    + "aca" + StringUtils.LS
    + "acg" + StringUtils.LS
    + "ccg" + StringUtils.LS
    ;

  /**
   * @throws IOException
   *
   */
  public void testIndel() throws IOException {
    checkSubstitute(0, 1, "acg" + StringUtils.LS, 0);
    checkSubstitute(1, 10, EXPECTED_INDEL, 0);
    checkSubstitute(2, 37, null, 0);
  }

  /**
   * @throws IOException
   *
   */
  public void testIndelAppend() throws IOException {
    checkSubstitute(0, 1, "acg" + StringUtils.LS, 1);
  }

}

