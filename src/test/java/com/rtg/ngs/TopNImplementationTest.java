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

import junit.framework.TestCase;

/**
 * Test class
 */
public class TopNImplementationTest extends TestCase {

  protected UptoNStore getTopNImplementation(final int numReads, final int numTemplateSeqs, final int n, final int templateMaxLength) {
    return new TopNImplementation(numReads, numTemplateSeqs, n, templateMaxLength, 50);
  }

  public void testToString() {
    Diagnostic.setLogStream();
    final UptoNStore topn = getTopNImplementation(4, 17, 5, 1000);
    assertEquals("TopNImplementation", topn.toString());
  }
  public void testDoesNotOverflowScoreIndelInTheResultCountArray() {
    Diagnostic.setLogStream();
    final UptoNStore topn = getTopNImplementation(5, 10, 50, 100);
    for (int i = 0; i < 255; ++i) {
      topn.process(1, false, 4, 1, 77);
    }
      topn.process(1, false, 4, 1, 5);
  }

  public void testOver255() {
    final UptoNStore topn = getTopNImplementation(5, 5, 1000, 2000);
    for (int i = 0; i < 1000; ++i) {
      topn.process(0, false, 2, 17, i % 50);
    }
    final MatchResult mr = new MatchResult(0);
    topn.setResults(mr, 2);
    assertEquals(1000, mr.size());
    for (int i = 0; i < mr.size(); ++i) {
      assertEquals(17, mr.getPosition(i));
    }
  }
  public void test() {
    Diagnostic.setLogStream();
    final UptoNStore topn = getTopNImplementation(4, 17, 5, 1000);
    topn.process(12, true, 0, 1, 5);
    topn.process(12, true, 0, 2, 7);
    topn.process(12, true, 0, 3, 2);
    topn.process(12, true, 0, 4, 2);
    topn.process(12, true, 0, 5, 7);
    topn.process(12, true, 0, 6, 5);
    topn.process(12, true, 0, 7, 9);
    topn.process(12, true, 0, 8, 1);
    topn.process(12, true, 0, 9, 0);
    topn.process(12, true, 0, 10, 10);
    topn.process(12, true, 0, 11, 6);
    topn.process(12, true, 0, 12, 5);

    topn.process(12, false, 1, 1, 6);
    topn.process(12, false, 1, 2, 4);
    topn.process(12, false, 1, 3, 16);
    topn.process(12, false, 1, 4, 7);
    topn.process(12, false, 1, 5, 17);
    topn.process(12, false, 1, 6, 10);
    topn.process(12, false, 1, 7, 3);

    topn.process(12, false, 2, 1, 3);
    topn.process(12, false, 2, 2, 3);
    topn.process(12, false, 2, 3, 3);
    topn.process(12, false, 2, 4, 3);
    topn.process(12, false, 2, 5, 3);
    topn.process(12, false, 2, 6, 3);
    topn.process(12, false, 2, 7, 3);
    topn.process(12, false, 2, 8, 3);

    topn.process(12, false, 3, 999, 4);
    topn.process(12, false, 3, 998, 2);

    MatchResult results = new MatchResult(2);
    topn.setResults(results, 0);
    assertEquals(4, results.size());
    checkReadResults(results, 0, 12, 9, true, 0);
    checkReadResults(results, 1, 12, 8, true, 0);
    checkReadResults(results, 2, 12, 3, true, 0);
    checkReadResults(results, 3, 12, 4, true, 0);

    results = new MatchResult(0);
    topn.setResults(results, 1);
    assertEquals(5, results.size());
    checkReadResults(results, 0, 12, 7, false, 1);
    checkReadResults(results, 1, 12, 2, false, 1);
    checkReadResults(results, 2, 12, 1, false, 1);
    checkReadResults(results, 3, 12, 4, false, 1);
    checkReadResults(results, 4, 12, 6, false, 1);

    results = new MatchResult(0);
    topn.setResults(results, 2);
    assertEquals(0, results.size());

    results = new MatchResult(0);
    topn.setResults(results, 3);
    assertEquals(2, results.size());
    checkReadResults(results, 0, 12, 998, false, 3);
    checkReadResults(results, 1, 12, 999, false, 3);
  }

  private static void checkReadResults(MatchResult result, int index, int templateId, int position, boolean reverse, int encReadId) {
    assertEquals(templateId, result.getTemplateId(index));
    assertEquals(position, result.getPosition(index));
    assertEquals(reverse, result.isReverse(index));
    assertEquals(encReadId, result.getEncodedReadId(index));
  }

}
