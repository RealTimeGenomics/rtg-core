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
package com.rtg.alignment;

import com.rtg.mode.DnaUtils;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class LoggingOnlyEditDistanceTest extends TestCase {

  /**
   * To test alignment speed, set repeatCount above to 100,000 then run this
   * test.
   */
  public void testEmpty() {

    MemoryPrintStream ms = new MemoryPrintStream();
    Diagnostic.setLogStream(ms.printStream());
    final byte[] read = new byte[0];
    final byte[] template = new byte[0];
    UnidirectionalEditDistance lo = new LoggingOnlyEditDistance();
    int[] res = lo.calculateEditDistance(read, 0, template, 0, 0, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    assertNull(res);
    // Logging only Test does just that?
    TestUtils.containsAll(ms.toString(), "LoggingOnlyEditDistance  1:", "crc=0", "read.length=0", "rlen=0", "zeroStart=0", "template.length=0",
        "maxScore=0", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

    lo = new LoggingOnlyEditDistance(3, "yo");
    ms = new MemoryPrintStream();
    Diagnostic.setLogStream(ms.printStream());
    res = lo.calculateEditDistance(DnaUtils.encodeString("AAAGCCCG"), 7, DnaUtils.encodeString("TTAAAGGCCGTT"), 1, 2, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    assertNull(res);
    TestUtils.containsAll(ms.toString(), "LoggingOnlyEditDistance yo 1:", "crc=3036913", "read.length=8", "rlen=7", "zeroStart=1", "template.length=12",
        "maxScore=2", "ED: AAAGCCCG TTAAAGGCCGTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN 20");

    ms = new MemoryPrintStream();
    Diagnostic.setLogStream(ms.printStream());
    res = lo.calculateEditDistance(DnaUtils.encodeString("AAAGCCCG"), 7, DnaUtils.encodeString("TTTTTGGGGGAAAAACCCCCTTAAAGGCCGTT"), 21, 2, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    assertNull(res);
    TestUtils.containsAll(ms.toString(), "ED: AAAGCCCG TTTTGGGGGAAAAACCCCCTTAAAGGCCGTTNNNNNNNNNNNNNNNN 20");
    Diagnostic.setLogStream();
  }



}
