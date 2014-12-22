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
    Graph g = GraphMapCliTest.makeGraph(2, new String[]{"AAACCT", "CTTTATATA"}, new long[][]{{1, 2}});
    MemoryPrintStream mps = new MemoryPrintStream();
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
