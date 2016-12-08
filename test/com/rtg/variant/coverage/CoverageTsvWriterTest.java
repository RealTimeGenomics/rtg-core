/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.coverage;

import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

public class CoverageTsvWriterTest extends TestCase {

  public void testTsvWriter() throws IOException {
    final MemoryPrintStream out = new MemoryPrintStream();
    final CoverageTsvWriter cw = new CoverageTsvWriter(out.outputStream());
    cw.init();
    final CoverageState cs = new CoverageState("hello", new byte[] {1, 2, 3, 4}, true);
    for (int i = 0; i < 4; ++i) {
      cs.incrementIH(i, 1);
      cs.incrementIH(i, 2);
    }
    cs.incrementIH(2, 1);
    cs.incrementIH(1, 2);
    for (int i = 1; i < 4; ++i) {
      cw.finalCoveragePosition("hello", i, cs.getIH1(i), cs.getIHgt1(i), cs.getScore(i));
    }
    final String actual = out.toString();
    TestUtils.containsAll(actual,
        CoverageBedWriter.VERSION_STRING,
        "#RUN-ID",
        "#Version",
        "Coverage output " + CoverageBedWriter.COVERAGE_OUTPUT_VERSION,
        "#sequence\tposition\tunique-count\tambiguous-count\tscore",
        "hello\t1\t1\t2\t2.00" + StringUtils.LS,
        "hello\t2\t2\t1\t2.50" + StringUtils.LS,
        "hello\t3\t1\t1\t1.50" + StringUtils.LS
    );
  }
}
