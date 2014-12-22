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

package com.rtg.util.arithcode;

import com.rtg.util.array.byteindex.ByteChunks;

import junit.framework.TestCase;

/**
 * Test for both <code>OutputBytes</code> and <code>InputBytes</code>
 */
public class BytesTest extends TestCase {

  private void write(OutputBytes ob, String str) {
    for (int i = 0; i < str.length(); i++) {
      ob.writeBit(str.charAt(i) == '1');
    }
  }

  private String read(InputBytes ib, int bits) {
    ib.integrity();
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < bits; i++) {
      sb.append(ib.readBit() ? "1" : "0");
    }
    assertFalse(ib.readBit());
    return sb.toString();
  }

  public void test() {
    final ByteChunks bc = new ByteChunks(0);
    final OutputBytes ob = new OutputBytes(bc);

    final String s1 = "00101";
    write(ob, s1);
    assertEquals(1, ob.endBlock());
    final String s2 = "00101111" + "1100011";
    write(ob, s2);
    assertEquals(3, ob.endBlock());
    ob.close();

    assertEquals(3, bc.length());
    assertEquals(0x28, bc.get(0));
    assertEquals(0x2F, bc.get(1));
    assertEquals(0xC6, bc.get(2));

    assertEquals(s1 + "000", read(new InputBytes(bc, 0, 1), 8));
    assertEquals(s2 + "0", read(new InputBytes(bc, 1, 3), 16));

  }

  //test empty and finish on byte boundary
  public void test1() {
    final ByteChunks bc = new ByteChunks(0);
    final OutputBytes ob = new OutputBytes(bc);

    write(ob, "");
    assertEquals(0, ob.endBlock());
    final String s2 = "00101111";
    write(ob, s2);
    assertEquals(1, ob.endBlock());
    ob.close();

    assertEquals(1, bc.length());
    assertEquals(0x2F, bc.get(0));

    assertEquals("0", read(new InputBytes(bc, 0, 0), 1));
    assertEquals(s2, read(new InputBytes(bc, 0, 1), 8));
  }

}
