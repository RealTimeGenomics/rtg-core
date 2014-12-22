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
public class SimpleOutputTest extends TestCase {

  public void test() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final SimpleOutput out = new SimpleOutput(bos, 1);
    assertEquals("SimpleOutput", out.toString());
    out.integrity();
    out.out("id", 1, 3, 2, 4);
    out.close();
    final String exp = "id\t0\t2\tcnv\t1.000\t1.414" + StringUtils.LS;
    assertEquals(exp, bos.toString());
  }

}
