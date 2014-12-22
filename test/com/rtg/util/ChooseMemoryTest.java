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
package com.rtg.util;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import junit.framework.TestCase;


/**
 */
public class ChooseMemoryTest extends TestCase {

  public void test() {
    final PrintStream o = System.out;
    final PrintStream e = System.err;
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final ByteArrayOutputStream err = new ByteArrayOutputStream();
    final PrintStream pout = new PrintStream(out);
    final PrintStream perr = new PrintStream(err);
    System.setOut(pout);
    System.setErr(perr);
    ChooseMemory.main(new String[] {});
    ChooseMemory.main(new String[] {"0"});
    ChooseMemory.main(new String[] {"101"});
    ChooseMemory.main(new String[] {"100"});
    ChooseMemory.main(new String[] {});
    pout.flush();
    perr.flush();
    final String outString = out.toString();
    final String errString = err.toString();
    assertTrue(outString.contains("m"));
    TestUtils.containsAll(errString, "Percentage must be greater than 0."
                                   , "Percentage must be less than or equal to 100.");
    pout.close();
    perr.close();
    System.setOut(o);
    System.setErr(e);
  }
}
