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
package com.rtg.index.hash.ngs.general;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import junit.framework.TestCase;


/**
 */
public class MaskAnalysisTest extends TestCase {

  public void test() throws Exception {
    final PrintStream oldOut = System.out;
    try {
      final ByteArrayOutputStream os = new ByteArrayOutputStream();
      final PrintStream newOut = new PrintStream(os);
      System.setOut(newOut);
      try {
        MaskAnalysis.main(new String[] {"6", "2", "1", "3"});
      } finally {
        newOut.close();
      }
      final String s = os.toString();
      //System.err.println(s);
      assertEquals("read=6substitutions=2indels=1indelLength=3" + LS
                   + "w\tbuild\tsearch\tsearch1\tsearch2\tmiss" + LS
                   + "1\t3\t#285714288.71#\t3\t#285714285.71#\t0.01%" + LS
                   + "2\t3\t#81632656.06#\t3\t#81632653.06#\t0.01%" + LS
                   + "3\t10\t#23323659.16#\t44\t#23323615.16#\t0.01%" + LS
                   + "4\t15\t6663968.05\t78\t6663890.05\t0.01%" + LS, s.replace(" ", ""));
    } finally {
      System.setOut(oldOut);
    }
  }

}
