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
package com.rtg.sam;

import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class SamFilterChainTest extends TestCase {

  public void testChaining() {
    final SamFilter filterAs = new SamFilter() {
      @Override
      public boolean acceptRecord(SAMRecord rec) {
        return !"a".equals(rec.getReadName());
      }
    };

    final SamFilter filterBs = new SamFilter() {
      @Override
      public boolean acceptRecord(SAMRecord rec) {
        return !"b".equals(rec.getReadName());
      }
    };

    final SamFilterChain chain = new SamFilterChain(filterAs, filterBs);

    final SAMRecord r1 = new SAMRecord(null);
    r1.setReadName("a");
    final SAMRecord r2 = new SAMRecord(null);
    r2.setReadName("b");
    final SAMRecord r3 = new SAMRecord(null);
    r3.setReadName("c");

    assertFalse(chain.acceptRecord(r1));
    assertFalse(chain.acceptRecord(r2));
    assertTrue(chain.acceptRecord(r3));
  }
}
