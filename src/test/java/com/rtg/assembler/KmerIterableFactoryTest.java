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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.StringUtils;
import com.rtg.util.iterators.Transform;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class KmerIterableFactoryTest extends TestCase {
  private static final Transform<Kmer, String> KMER2STRING = new Kmer2String();
  private static class Kmer2String extends Transform<Kmer, String> {
    @Override
    public String trans(Kmer x) {
      return x.toString();
    }
  }

  private static <X> void check(Iterator<X> exp, Iterator<X> actual) {
    while (exp.hasNext()) {
      assertTrue(actual.hasNext());
      final X e = exp.next();
      final X a = actual.next();
      //System.err.println(e + ":" + a);
      assertEquals(e, a);
    }
    assertFalse(actual.hasNext());
  }

  public void test() throws IOException {
    final File tmpDir =  ReaderTestUtils.getDNADir(">read1" + StringUtils.LS + "accgttc" + StringUtils.LS);
    try {
      try (ReadPairSource source = ReadPairSource.makeSource(tmpDir, LongRange.NONE)) {
        final KmerIterableFactoryInterface factory = new KmerIterableFactory(Collections.singletonList(source), StringKmer.factory(), 4);
        try (KmerIterable kmers = factory.makeIterable()) {
          final Iterator<String> it = KMER2STRING.trans(kmers.iterator());
          final Iterator<String> exp = Transform.array2Iterator(new String[]{"ACCG", "CCGT", "CGTT", "GTTC"});
          check(exp, it);
        }
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
