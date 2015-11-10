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
  protected void tearDown() throws Exception {
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
    for (int i = 0; i < 3; i++) {
      assertEquals(0, ((MockIndex) is.get(i)).getTimesFrozen());
    }
    final MemoryPrintStream baos = new MemoryPrintStream();
    Diagnostic.setLogStream(baos.printStream());
    is.freeze(2);
    for (int i = 0; i < 3; i++) {
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
