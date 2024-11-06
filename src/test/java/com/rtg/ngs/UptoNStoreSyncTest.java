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
package com.rtg.ngs;

import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class UptoNStoreSyncTest extends TopNImplementationTest {

  @Override
  protected UptoNStore getTopNImplementation(final int numReads, final int numTemplateSeqs, final int n, final int templateMaxLength) {
    return new UptoNStoreSync(new TopNImplementation(numReads, numTemplateSeqs, n, templateMaxLength, 50));
  }

  @Override
  public void testToString() {
    Diagnostic.setLogStream();
    final UptoNStore topn = getTopNImplementation(4, 17, 5, 1000);
    assertEquals("TopNImplementationSync", topn.toString());
  }

  public void testSyncId() {
    assertEquals(0, UptoNStoreSync.syncId(0));
    assertEquals(0, UptoNStoreSync.syncId(65536));
    assertEquals(1, UptoNStoreSync.syncId(1));
    assertEquals(1, UptoNStoreSync.syncId(65536 + 1));
  }
}
