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
package com.rtg.util.io.bzip2;

import java.util.zip.CRC32;

/**
 * Test class
 */
public class CRCTest extends CBZip2InputStreamTest {

  public void testCRC() {
    //it seems that the bzip2 CRC algorithm uses a reversed bit order to that of the standard Crc32 used by gzip
    final CRC crc = new CRC();
    final CRC32 exp = new CRC32();
    crc.updateCRC(42);
    crc.updateCRC(42);
    crc.updateCRC(45);
    exp.update(84); //byte 42 reversed
    exp.update(84);
    exp.update(180); //byte 45 reversed

    assertEquals(Long.toBinaryString(exp.getValue()), new StringBuilder(Long.toBinaryString(crc.getFinalCRC() & 0xffffffffL)).reverse().toString());

    final CRC crc2 = new CRC();
    crc2.updateCRC(42, 2);
    crc2.updateCRC(45, 1);

    assertEquals(Long.toBinaryString(exp.getValue()), new StringBuilder(Long.toBinaryString(crc2.getFinalCRC() & 0xffffffffL)).reverse().toString());
  }

}
