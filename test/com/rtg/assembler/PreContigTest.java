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

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 *         Date: 10/05/12
 *         Time: 4:31 PM
 */
public class PreContigTest extends TestCase {
  public void testPreContig() {
    PreContig contig = new PreContig(3, new ByteKmer("GTC"), 7);
    assertEquals(7, contig.mKmerCount);
    assertEquals("GTC", contig.mContig.toString());
    contig.extend(true, new ByteKmer("TCT"), 2);
    assertEquals(9, contig.mKmerCount);
    assertEquals("GTCT", contig.mContig.toString());
    contig.extend(true, new ByteKmer("CTC"), 1);
    assertEquals(10, contig.mKmerCount);
    assertEquals("GTCTC", contig.mContig.toString());
    contig.extend(false, new ByteKmer("AGT"), 7);
    assertEquals(17, contig.mKmerCount);
    assertEquals("AGTCTC", contig.mContig.toString());

    assertEquals("PreContig: id=3 sequence:AGTCTC" + StringUtils.LS, contig.toString());
  }

}
