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
