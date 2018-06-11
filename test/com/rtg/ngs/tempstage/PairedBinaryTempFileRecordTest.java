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
package com.rtg.ngs.tempstage;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class PairedBinaryTempFileRecordTest extends TestCase {
  public void testBasics() throws IOException {
    final BinaryTempFileRecord bar = new BinaryTempFileRecord(true, false, false, false);
    bar.setCigarString(new byte[0]);
    bar.setMdString(new byte[0]);
    bar.setComboScore(5);
    bar.setMatePosition(21);
    bar.setTemplateLength(100);
    bar.setSamFlags((byte) 65);
    check(bar);
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (TempRecordWriterNio writer = new TempRecordWriterNio(baos)) {
      writer.writeRecord(bar);
      bar.setSentinelRecord();
      writer.writeRecord(bar);
    }
    try (TempRecordReaderNio dis = new TempRecordReaderNio(new ByteArrayInputStream(baos.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, false))) {
      final BinaryTempFileRecord rec;
      assertNotNull(rec = dis.readRecord());
      check(rec);
      assertNull(dis.readRecord());
    }
  }

  private void check(BinaryTempFileRecord bar) {
    assertEquals(100, bar.getTemplateLength());
    assertEquals(5,  bar.getComboScore());
    assertEquals(21, bar.getMatePosition());
    assertEquals(65, bar.getSamFlags());
    assertTrue(bar.isReadPaired());
    assertFalse(bar.isMateReverseStrand());
    assertTrue(bar.isFirstOfPair());
    assertFalse(bar.isSecondOfPair());
  }
}
