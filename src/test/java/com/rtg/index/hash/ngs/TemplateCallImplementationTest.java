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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.index.Finder;
import com.rtg.index.Index;
import com.rtg.index.IndexSet;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest;
import com.rtg.index.hash.ngs.instances.SplitL4w2s1e1;
import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.DefaultOutputProcessor;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsTestUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class TemplateCallImplementationTest extends TestCase {

  @Override
  protected void tearDown() {
    Diagnostic.setLogStream();
  }

  private NgsHashFunction getHashFunction(final IndexSet indexes, final TemplateCall tc) {
    final SplitL4w2s1e1 hf = new SplitL4w2s1e1(new ReadCallImplementation(indexes), tc);
    hf.integrity();
    AbstractSplitTest.encode(hf, "acta");
    return hf;
  }

  public final void test() throws IOException {
    Diagnostic.setLogStream();
    final ByteArrayOutputStream sb = new ByteArrayOutputStream();
    final Index[] indexes = new Index[3];
    for (int i = 0; i < indexes.length; ++i) {
      indexes[i] = new IndexMock(sb, i) {
          @Override
          public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
            finder.found(mIndex);
          }
        };
    }

    final NgsFilterParams filterParams = NgsFilterParams.builder().errorLimit(3).create();
    final NgsTestUtils.OverriddenNgsOutputParams outParams = new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(sb).filterParams(filterParams));
    final File file = ReaderTestUtils.getDNADir(">t\na");
    final SequenceParams buildParams = SequenceParams.builder().directory(file).create();
    final NgsParams params = NgsParams.builder().outputParams(outParams).buildFirstParams(buildParams).create();
    final DefaultOutputProcessor outPr = new DefaultOutputProcessor(params);
    final IndexSet indexSet = new IndexSet(indexes);
    final TemplateCallImplementation tci = new TemplateCallImplementation(params, 4, indexSet, outPr);
    final NgsHashFunction hashFunction = getHashFunction(indexSet, tci);
    hashFunction.setReadSequences(8);
    hashFunction.setValues(0, false);

    tci.setReverse(false);
    tci.set(123, 9);
    tci.templateCall(4, 3, 0);
    tci.done();

    tci.setReverse(false);
    tci.set(123, 10);
    tci.templateCall(5, 3, 0);
    tci.templateCall(5, 4, 0);
    tci.templateCall(5, 4, 1);
    tci.done();

    tci.setReverse(true);
    tci.set(234, 6);
    tci.templateCall(3, 3, 0);
    tci.templateCall(3, 4, 0);
    tci.templateCall(3, 4, 1);
    tci.done();

    outPr.finish();
    //System.err.println(sb.toString());
    TestUtils.containsAll(sb.toString(), "123\tF\t0\t2", "123\tF\t0\t3", "123\tF\t1\t3", "234\tR\t0\t1");
    assertEquals("TemplateCallImplementation", tci.toString());
    assertTrue(FileHelper.deleteAll(file));
  }

  public final void testClone() throws IOException, CloneNotSupportedException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(baos);
    try {
      Diagnostic.setLogStream(ps);
      final ByteArrayOutputStream sb = new ByteArrayOutputStream();
      final Index[] indexes = new Index[3];
      for (int i = 0; i < indexes.length; ++i) {
        indexes[i] = new IndexMock(sb, i) {
            @Override
            public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
              finder.found(mIndex);
            }
          };
      }

      final NgsFilterParams filterParams = NgsFilterParams.builder().errorLimit(3).create();
      final NgsOutputParamsBuilder outBuilder = NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(sb).filterParams(filterParams);
      final NgsTestUtils.OverriddenNgsOutputParams outParams = new NgsTestUtils.OverriddenNgsOutputParams(outBuilder);
      final File file = ReaderTestUtils.getDNADir(">t\na");
      final SequenceParams buildParams = SequenceParams.builder().directory(file).create();
      final NgsParams params = NgsParams.builder().outputParams(outParams).buildFirstParams(buildParams).create();
      final DefaultOutputProcessor outPr = new DefaultOutputProcessor(params);
      final IndexSet indexSet = new IndexSet(indexes);
      final TemplateCallImplementation tca = new TemplateCallImplementation(params, 4, indexSet, outPr);
      assertTrue(tca.integrity());
      final TemplateCallImplementation tci = tca.clone();
      final NgsHashFunction hashFunction = getHashFunction(indexSet, tci);
      assertTrue(tci.integrity());
      hashFunction.setReadSequences(8);
      hashFunction.setValues(0, false);

      tci.setReverse(false);
      tci.set(123, 9);
      tci.templateCall(4, 3, 0);
      tci.done();

      tci.setReverse(false);
      tci.set(123, 10);
      tci.templateCall(5, 3, 0);
      tci.templateCall(5, 4, 0);
      tci.templateCall(5, 4, 1);
      tci.done();

      tci.setReverse(true);
      tci.set(234, 6);
      tci.templateCall(3, 3, 0);
      tci.templateCall(3, 4, 0);
      tci.templateCall(3, 4, 1);
      tci.done();

      outPr.finish();
      TestUtils.containsAll(sb.toString(), "123\tF\t0\t2", "123\tF\t0\t3", "123\tF\t1\t3", "234\tR\t0\t1");
      assertTrue(FileHelper.deleteAll(file));
      tci.logStatistics();
    } finally {
      ps.flush();
      ps.close();
    }
    assertTrue(baos.toString().contains("TemplateCall Hits 7"));
    assertTrue(baos.toString().contains("TemplateCall Processed 4"));
  }
}

