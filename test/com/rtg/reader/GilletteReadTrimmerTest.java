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

package com.rtg.reader;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 */
public class GilletteReadTrimmerTest extends TestCase {

  private void check(String raw, int expLenq15w15, int expLenq15w30) {
    final GilletteReadTrimmer grt15 = new GilletteReadTrimmer(15, 15);
    final GilletteReadTrimmer grt30 = new GilletteReadTrimmer(30, 15);

    final byte[] quals = DnaUtils.fastqToPhred(raw);

    assertEquals(expLenq15w15, grt15.getTrimPosition(quals, raw.length()));
    assertEquals(expLenq15w30, grt30.getTrimPosition(quals, raw.length()));

  }

  public void testExamples() {
    check(";88/;5=8??<??;?=;2222*274=???0?;????????3????3????>8?7>>>=>?>?>==:97*744*7799/=6=<975)'''+,',/", 94, 94);
    check("??<?89.43---+3,79=?5?;??=??????????>??;>>?9;8441,,+14479;9976,2//-,,,-/),13++.43/12", 73, 83);
    check("::::18:5:4:458:::':::8664/.(1,,''.22(262/.", 35, 42);
    check("3.1-'''''))'11+2.6811::88:6::484,440111*))'--''**1+4.55.32(*--''''-')*1+0+23++45661885", 0, 56);
    check(">???>??>???", 0, 0);  //shorter than the window length!
  }

}

