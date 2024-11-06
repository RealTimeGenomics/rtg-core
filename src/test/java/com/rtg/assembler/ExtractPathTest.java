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

package com.rtg.assembler;

import com.rtg.assembler.graph.Graph;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class ExtractPathTest extends AbstractCliTest {
  public void testHelp() {
    checkHelp("rtg ExtractPath"
        , "input graph directory"
        , "-p,", "--path", "path that will be traced through the graph and output"
        , "-k,", "--kmer-size", "size of the kmer the graph was constructed with"
        , "-s,", "--start", "trim the first sequence up to this position"
        , "-e,", "--end", "trim the final sequence at this position"
    );
  }

  public void testExtract() {
    final Graph g = GraphMapCliTest.makeGraph(2, new String[]{"AAACCT", "CTTTATATA"}, new long[][]{{1, 2}});
    final MemoryPrintStream mps = new MemoryPrintStream();
    ExtractPath.run(mps.outputStream(), g, 1, 3, 2, 6);
    assertEquals("ACCT|TTATA" + StringUtils.LS, mps.toString());

    mps.reset();
    ExtractPath.run(mps.outputStream(), g, 1, 3, -1, -1);
    assertEquals("AAACCT|TTATATA" + StringUtils.LS, mps.toString());

    mps.reset();

    ExtractPath.run(mps.outputStream(), g, 1, 3, 0, -1);
    assertEquals("AAACCT|TTATATA" + StringUtils.LS, mps.toString());
  }

  @Override
  protected AbstractCli getCli() {
    return new ExtractPath();
  }
}
