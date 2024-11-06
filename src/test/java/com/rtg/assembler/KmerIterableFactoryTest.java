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
