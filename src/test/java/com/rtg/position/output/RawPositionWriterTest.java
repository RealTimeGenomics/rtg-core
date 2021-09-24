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
package com.rtg.position.output;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.StringWriter;

import com.rtg.mode.UnidirectionalFrame;

/**
 * Tests the corresponding class.
 *
 */
public class RawPositionWriterTest extends AbstractPositionWriterTest {

  public void test() throws IOException {
    final StringWriter out = new StringWriter();
    final StringWriter ambiguous = new StringWriter();
    final StringWriter unmapped = new StringWriter();
    final PositionWriter w = new RawPositionWriter(out, unmapped);
    assertEquals(0, w.score(), 1E-12);
    w.endQuery(1);
    w.write(getRegion(2, 1, 32, 32.0), 43, UnidirectionalFrame.FORWARD, 44, 0 /*unused*/);
    w.endQuery(43);
    w.write(getRegion(2, 1, 32, 32.0), 53, UnidirectionalFrame.FORWARD, 44, 0);
    w.write(getRegion(2, 1, 32, 32.0), 53, UnidirectionalFrame.FORWARD, 44, 0);
    w.endQuery(53);
    w.endQuery(63);
    assertEquals("1" + LS + "63" + LS, unmapped.toString());
    final String exp = ""
      + "#query-id\tquery-frame\tquery-start\tquery-end\tsubject-id\tsubject-frame\tsubject-start\tsubject-end" + LS
      + "43\t\t33\t64\t2\t\t2\t33" + LS
      + "53\t\t33\t64\t2\t\t2\t33" + LS
      + "53\t\t33\t64\t2\t\t2\t33" + LS
      ;
    assertEquals(exp, out.toString());
    assertEquals("", ambiguous.toString());
    assertEquals(96.0, w.score());
  }
}

