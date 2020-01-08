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
package com.rtg.protein;

import java.io.IOException;

import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;

/**
 *
 */
public class TopEqualProteinImplementationTest extends TopEqualProteinOutputProcessorTest {

  /**
   */
  public TopEqualProteinImplementationTest(String name) {
    super(name);
  }

  public void testResults() throws IOException, InvalidParamsException {
    final TopEqualProteinImplementation impl = new TopEqualProteinImplementation(2, 3);
    impl.insertResult(getResult(0));
    impl.insertResult(getResult(1));
    impl.insertResult(getResult(1));
    impl.insertResult(getResult(2));
    impl.insertResult(getResult(2));
    impl.insertResult(getResult(2));
    impl.contains(0, 0, 1, 1) ;

    assertEquals(1, impl.resultCount(0));
    assertEquals(2, impl.resultCount(1));
    assertEquals(3, impl.resultCount(2));
  }

  public void testResults2() throws IOException, InvalidParamsException {
    final TopEqualProteinImplementation impl = new TopEqualProteinImplementation(2, 3);
    impl.insertResult(getResult(0));
    impl.insertResult(getResult(1));
    impl.insertResult(getResult(1));
    impl.insertResult(getResult(1));
    impl.contains(0, 0, 1, 1) ;
    assertEquals(1, impl.resultCount(0));
    assertEquals(3, impl.resultCount(1));
  }

  private static ProteinAlignmentResult getResult(int readId) throws InvalidParamsException, IOException {
    final SharedProteinResources resx = new SharedProteinResources(new ProteinScoringMatrix(), ReaderTestUtils.getReaderProteinMemory(">A\n"), ReaderTestUtils.getReaderDnaMemory(">A\n"), false);
    return new ProteinAlignmentResult(resx, 0, readId * 6, null, 0, true);
  }
}
