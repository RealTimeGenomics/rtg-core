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
