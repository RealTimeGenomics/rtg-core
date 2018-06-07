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

