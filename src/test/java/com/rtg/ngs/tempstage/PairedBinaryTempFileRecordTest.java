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
