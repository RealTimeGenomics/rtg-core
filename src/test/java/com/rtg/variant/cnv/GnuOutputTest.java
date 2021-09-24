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
package com.rtg.variant.cnv;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class GnuOutputTest extends TestCase {

  public void test() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final GnuOutput out = new GnuOutput(bos);
    assertEquals("GnuOutput", out.toString());
    out.out("id", 1, 5, 2, 4);
    out.out("id", 5, 6, 2, 4);
    out.close();
    final String exp = ""
      + "-0.25 0.500" + StringUtils.LS
      + "3.25 0.500" + StringUtils.LS
      + StringUtils.LS
      + "-0.25 -0.500" + StringUtils.LS
      + "-0.25 1.500" + StringUtils.LS
      + StringUtils.LS
      + "3.25 -0.500" + StringUtils.LS
      + "3.25 1.500" + StringUtils.LS
      + StringUtils.LS
      + "3.75 2.000" + StringUtils.LS
      + "4.25 2.000" + StringUtils.LS
      + StringUtils.LS
      + "3.75 2.000" + StringUtils.LS
      + "3.75 2.000" + StringUtils.LS
      + StringUtils.LS
      + "4.25 2.000" + StringUtils.LS
      + "4.25 2.000" + StringUtils.LS
      + StringUtils.LS
      ;
    assertEquals(exp, bos.toString());
  }

}
