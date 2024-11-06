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
package com.rtg.position.output;

import java.io.File;
import java.io.IOException;
import java.io.StringWriter;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Test the enumeration OutputType.
 */
public class OutputFormatTypeTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void test() {
    assertEquals("[SEGMENT, NGS]", OutputFormatType.values().toString());
  }


  private PositionParams makeParams(final OutputFormatType type) throws IOException {
    final ReaderParams srp = new MockReaderParams(0, 0, SequenceMode.UNIDIRECTIONAL.codeType());
    final ISequenceParams subjectParams = new MockSequenceParams(srp, SequenceMode.UNIDIRECTIONAL, 0, 0);
    final BuildParams buildParams = BuildParams.builder().windowSize(5).stepSize(2).sequences(subjectParams).create();

    final ReaderParams qrp = new MockReaderParams(0, 0, SequenceMode.BIDIRECTIONAL.codeType());
    final ISequenceParams queryParams = new MockSequenceParams(qrp, SequenceMode.BIDIRECTIONAL, 0, 0);
    final BuildParams searchParams = BuildParams.builder().windowSize(5).stepSize(1).sequences(queryParams).create();
    final PositionDistributionParams distr = new PositionDistributionParams(0.001, 0.009, 20, 0);
    final PositionOutputParams outParams = new PositionOutputParams(new File(""), type, distr, null, false, 10);
    return PositionParams.builder().hashCountThreshold(1000).buildParams(buildParams).searchParams(searchParams).outputParams(outParams).create();
  }

  public void testOutput() throws IOException {
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    for (final OutputFormatType type : OutputFormatType.values()) {
      //System.err.println(type);
      final PositionParams params = makeParams(type);
      final GapBucketsInfo bucketInfo = new GapBucketsInfo(params, 1);
      final PositionOutput output = type.output(params, bucketInfo, null, out, unmappedOut, null, params.hashCountThreshold());
      assertNotNull(output);
    }
  }

}
