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
    cw.finalCoveragePosition("hello", 1, 1, 2, 2.0);
    cw.finalCoveragePosition("hello", 2, 2, 1, 2.5);
    cw.finalCoveragePosition("hello", 3, 1, 1, 1.5);
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
