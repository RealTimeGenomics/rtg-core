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

import junit.framework.TestCase;

/**
 */
public class DeduplicatingNStoreTest extends TestCase {
  public void test() {
    final DeduplicatingNStore store = new DeduplicatingNStore(10, 5, 3, 102, 50);
    store.process(0, true, 0, 5, 0);
    store.process(0, true, 0, 400, 0);
    store.process(0, true, 0, 400, 0);
    store.process(0, true, 0, 5, 0);
    store.process(0, true, 0, 6, 0);
    store.process(0, true, 1, 5, 0);
    store.process(0, true, 1, 400, 0);
    store.process(0, true, 1, 400, 0);
    store.process(0, true, 1, 5, 0);
    store.process(0, true, 1, 7, 0);
    store.process(0, true, 1, 6, 0);
    final MatchResult r = new MatchResult(2);
    store.setResults(r, 0);
    store.setResults(r, 1);
    assertEquals(3, r.size());
    assertEquals(5, r.getPosition(0));
    assertEquals(true, r.isReverse(0));
    assertEquals(0, r.getEncodedReadId(0));
    assertEquals(0, r.getTemplateId(0));
    assertEquals(6, r.getPosition(1));
    assertEquals(true, r.isReverse(1));
    assertEquals(0, r.getEncodedReadId(1));
    assertEquals(0, r.getTemplateId(1));
    assertEquals(400, r.getPosition(2));
    assertEquals(true, r.isReverse(2));
    assertEquals(0, r.getEncodedReadId(2));
    assertEquals(0, r.getTemplateId(2));
  }
}
