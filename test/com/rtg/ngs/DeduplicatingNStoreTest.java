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
package com.rtg.ngs;

import junit.framework.TestCase;

/**
 */
public class DeduplicatingNStoreTest extends TestCase {
  public void test() {
    DeduplicatingNStore store = new DeduplicatingNStore(10, 5, 3, 102, 50);
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
    MatchResult r = new MatchResult(2);
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
