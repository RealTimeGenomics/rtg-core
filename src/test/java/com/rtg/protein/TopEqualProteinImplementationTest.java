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
package com.rtg.protein;

import java.io.IOException;

import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;

/**
 *
 */
public class TopEqualProteinImplementationTest extends TopEqualProteinOutputProcessorTest {

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
