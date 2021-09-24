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

import java.io.IOException;
import java.io.StringWriter;

import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallMock;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.TemplateCallMock;
import com.rtg.index.hash.ngs.protein.ProteinMask;

import junit.framework.TestCase;

/**
 * Tests the corresponding class
 */
public class ProteinTemplateHashLoopTest extends TestCase {

  public void testHashLoop1a() throws IOException {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final Skeleton sk = new Skeleton(4, 4, 0, 0, 1);
    final HashFunctionFactory f = ProteinMask.factory(sk);
    final NgsHashFunction hf = f.create(rcall, call);

    final ProteinTemplateHashLoop hashLoop1 = new ProteinTemplateHashLoop(5, 1, (ProteinMask) hf);
    hashLoop1.hashCall(1, 1);
  }

  public void testUnsupported() {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final Skeleton sk = new Skeleton(4, 4, 0, 0, 1);
    final HashFunctionFactory f = ProteinMask.factory(sk);
    final NgsHashFunction hf = f.create(rcall, call);

    final ProteinTemplateHashLoop hashLoop1 = new ProteinTemplateHashLoop(5, 1, (ProteinMask) hf);
    try {
      hashLoop1.hashCallBidirectional(1, 1, 1, 1);
      fail();
    } catch (UnsupportedOperationException e) {
      // Succeeded
      assertEquals("Not supported.", e.getMessage());
    }
  }


}
