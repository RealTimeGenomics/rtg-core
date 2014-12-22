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
package com.rtg.index.hash.ngs;


import java.io.ByteArrayOutputStream;

import com.rtg.index.Index;
import com.rtg.index.IndexSet;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class ReadCallImplementationTest extends TestCase {

  /**
   * Test method for {@link com.rtg.index.hash.ngs.ReadCallImplementation#readCall(int, long, int)}.
   */
  public final void testReadCall() {
    final ByteArrayOutputStream sb = new ByteArrayOutputStream();
    final Index[] indexes = new Index[3];
    for (int i = 0; i < indexes.length; i++) {
      indexes[i] = new IndexMock(sb, i);
    }
    final ReadCallImplementation rci = new ReadCallImplementation(new IndexSet(indexes));
    rci.readCall(1, 101L, 0);
    rci.readCall(2, 102L, 1);
    rci.readCall(3, 103L, 2);
    assertEquals(""
        + "add hash=101 index=0 id=1" + StringUtils.LS
        + "add hash=102 index=1 id=2" + StringUtils.LS
        + "add hash=103 index=2 id=3" + StringUtils.LS
        , sb.toString()
    );
  }

}

