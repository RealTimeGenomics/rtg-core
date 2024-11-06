/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
