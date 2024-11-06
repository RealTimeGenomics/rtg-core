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
package com.rtg.index;

import java.io.IOException;

import com.rtg.index.params.CreateParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.position.MockIndex;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class IndexSetTest extends TestCase {

  @Override
  protected void tearDown() {
    Diagnostic.setLogStream();
  }

  public void testSimpleCase() throws IOException {
    final Index[] indexes = {
        new MockIndex()
        , new MockIndex()
        , new MockIndex()
    };
    final IndexSet is = new IndexSet(indexes);
    assertEquals(3, is.size());
    for (int i = 0; i < 3; ++i) {
      assertEquals(0, ((MockIndex) is.get(i)).getTimesFrozen());
    }
    final MemoryPrintStream baos = new MemoryPrintStream();
    Diagnostic.setLogStream(baos.printStream());
    is.freeze(2);
    for (int i = 0; i < 3; ++i) {
      assertEquals(1, ((MockIndex) is.get(i)).getTimesFrozen());
    }
    final String str = baos.toString();
    //System.err.println(str);
    TestUtils.containsAll(str
        , "Start freeze job 0"
        , "Start freeze job 1"
        , "Start freeze job 2"
        , "Worker Thread Created - BuildFreeze-0 - 1/2"
        );

  }
  public void testComplexConstructor() throws IOException {
    final NgsParamsBuilder paramsBuilder = new NgsParamsBuilder().numberThreads(2).outputParams(NgsOutputParams.builder().create());
    final NgsParams params = paramsBuilder.create();
    final CreateParams.CreateParamsBuilder indexBuilder = new CreateParams.CreateParamsBuilder();
    final CreateParams indexParams = indexBuilder.create();

    final MemoryPrintStream baos = new MemoryPrintStream();
    Diagnostic.setLogStream(baos.printStream());
    final IndexSet is = new IndexSet(params, indexParams, 2);
    final int expectedLength = 2;
    assertEquals(expectedLength, is.size());
    assertNotNull(is.get(expectedLength - 1));
    TestUtils.containsAll(baos.toString()
        , "Start create job 0"
        , "Start create job " + (expectedLength - 1)
        , "maximum 2 threads"
        , "Worker Thread Created - CreateIndex-0"
        );
  }
}
