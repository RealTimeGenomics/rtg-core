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
  protected void tearDown() throws Exception {
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
    for (int i = 0; i < indexes.length; i++) {
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
    final NgsParams params = NgsParams.builder().hashCountThreshold(1000).outputParams(outParams).buildFirstParams(buildParams).create();
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
      for (int i = 0; i < indexes.length; i++) {
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
      final NgsParams params = NgsParams.builder().hashCountThreshold(1000).outputParams(outParams).buildFirstParams(buildParams).create();
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

