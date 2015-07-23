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
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.index.Finder;
import com.rtg.index.Index;
import com.rtg.index.IndexSet;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest;
import com.rtg.index.hash.ngs.instances.SplitL4w2s1e1;

import junit.framework.TestCase;

/**
 */
public class PcrTemplateCallImplementationTest extends TestCase {

  private static final int IX = 3;

  private NgsHashFunction getHashFunction(final IndexSet indexes, final TemplateCall tc) {
    final SplitL4w2s1e1 hf = new SplitL4w2s1e1(new ReadCallImplementation(indexes), tc);
    hf.integrity();
    AbstractSplitTest.encode(hf, "acta");
    return hf;
  }

  public final void testClone() throws IOException, CloneNotSupportedException {
    final OutputStream sb = new ByteArrayOutputStream();
    final Index[] indexes = new Index[IX];
    for (int i = 0; i < indexes.length; i++) {
      indexes[i] = new IndexMock(sb, i) {
          @Override
          public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
          finder.found(mIndex);
        }
      };
    }

    final byte[] dist = new byte[8];
    for (int i = 0; i < dist.length; i++) {
      dist[i] = 4;
    }
    final IndexSet indexSet = new IndexSet(indexes);
    final PcrTemplateCallImplementation tca = new PcrTemplateCallImplementation(indexSet, dist);
    assertTrue(tca.integrity());
    final PcrTemplateCallImplementation tci = tca.clone();
    final NgsHashFunction hashFunction = getHashFunction(indexSet, tci);
    assertTrue(tci.integrity());
    hashFunction.setReadSequences(8);
    tci.setHashFunction(hashFunction);

    tci.templateCall(0, 3, 0);
    tci.done();

    tci.templateCall(0, 3, 1);
    tci.templateCall(0, 4, 2);
    tci.done();
    //System.err.println(Arrays.toString(dist));
    for (int i = 0; i < IX; i++) {
      assertEquals(2, dist[i]);
    }
    for (int i = IX; i < 8; i++) {
      assertEquals(4, dist[i]);
    }
  }

  public final void test() throws IOException {
    final OutputStream sb = new ByteArrayOutputStream();
    final Index[] indexes = new Index[IX];
    for (int i = 0; i < indexes.length; i++) {
      indexes[i] = new IndexMock(sb, i) {
          @Override
          public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
          finder.found(mIndex);
        }
      };
    }

    final byte[] dist = new byte[8];
    for (int i = 0; i < dist.length; i++) {
      dist[i] = 4;
    }
    final IndexSet indexSet = new IndexSet(indexes);
    final PcrTemplateCallImplementation tci = new PcrTemplateCallImplementation(indexSet, dist);
    final NgsHashFunction hashFunction = getHashFunction(indexSet, tci);
    hashFunction.setReadSequences(8);
    tci.setHashFunction(hashFunction);

    tci.templateCall(0, 3, 0);
    tci.done();

    tci.templateCall(0, 3, 1);
    tci.templateCall(0, 4, 2);
    tci.done();
    //System.err.println(Arrays.toString(dist));
    for (int i = 0; i < IX; i++) {
      assertEquals(2, dist[i]);
    }
    for (int i = IX; i < 8; i++) {
      assertEquals(4, dist[i]);
    }
    assertEquals("PcrTemplateCallImplementation", tci.toString());
  }
}

