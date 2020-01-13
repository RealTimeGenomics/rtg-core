/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
