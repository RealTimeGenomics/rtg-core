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

package com.rtg.reader;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayInputStream;

import org.junit.Test;

/**
 */
public class FastqPairIteratorTest {
  @Test
  public void empty() {
    final FastqPairIterator fastqPairIterator = new FastqPairIterator(emptyFastqIterator(), emptyFastqIterator());
    assertFalse(fastqPairIterator.hasNext());

  }

  private static final String LEFT_FASTQ = "@name\nAAAAAAAA\n+name\nBBBBBBBB\n"
    + "@second\nGGGG\n+second\nAAAA\n";
  private static final String RIGHT_FASTQ = "@name\nCCCCCCCC\n+name\nBBBBBBBB\n"
    + "@second\nTTTT\n+second\nAAAA\n";

  @Test
  public void full() {
    final FastqSequenceDataSource leftSource = new FastqSequenceDataSource(new ByteArrayInputStream(LEFT_FASTQ.getBytes()), QualityFormat.SOLEXA);
    final FastqSequenceDataSource rightSource = new FastqSequenceDataSource(new ByteArrayInputStream(RIGHT_FASTQ.getBytes()), QualityFormat.SOLEXA);
      final FastqPairIterator fastqIterator = new FastqPairIterator(new FastqIterator(leftSource), new FastqIterator(rightSource));
      assertTrue(fastqIterator.hasNext());
      FastqPair pair = fastqIterator.next();
      FastqIteratorTest.checkRead(pair.r1(), "name", "AAAAAAAA");
      FastqIteratorTest.checkRead(pair.r2(), "name", "CCCCCCCC");
      assertTrue(fastqIterator.hasNext());
      pair = fastqIterator.next();
      FastqIteratorTest.checkRead(pair.r1(), "second", "GGGG");
      FastqIteratorTest.checkRead(pair.r2(), "second", "TTTT");
      assertFalse(fastqIterator.hasNext());
    }


  FastqIterator emptyFastqIterator() {

    final FastqSequenceDataSource fastqSequenceDataSource = new FastqSequenceDataSource(new ByteArrayInputStream(new byte[0]), QualityFormat.SOLEXA);
    return new FastqIterator(fastqSequenceDataSource);
  }

}
