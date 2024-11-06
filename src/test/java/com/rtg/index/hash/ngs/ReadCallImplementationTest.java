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
    for (int i = 0; i < indexes.length; ++i) {
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

