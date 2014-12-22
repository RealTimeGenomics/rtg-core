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

package com.rtg.tabix;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 * Test class
 */
public class BrLineReaderTest extends TestCase {

  private static final String EXPECTED = ""
      + "simulatedSequence19\t583\t.\tA\tT\t.\tPASS\t.\tGT:DP:RE:GQ:RS\t0/1:23:0.458:185.0:A,11,0.219,T,12,0.239" + StringUtils.LS
      + "simulatedSequence19\t637\t.\tG\tC\t.\tPASS\t.\tGT:DP:RE:GQ:RS\t1/0:27:0.537:53.0:C,7,0.139,G,20,0.398" + StringUtils.LS
      + "simulatedSequence19\t737\t.\tG\tC\t.\tPASS\t.\tGT:DP:RE:GQ:RS\t1/1:27:0.537:74.0:C,26,0.517,T,1,0.020" + StringUtils.LS;

  public void testSomeMethod() throws IOException {
    try (BrLineReader br = new BrLineReader(new BufferedReader(new StringReader(EXPECTED)))) {
      int count = 0;
      while (br.readLine() != null) {
        count++;
      }
      assertEquals(3, count);
    }
  }

}
