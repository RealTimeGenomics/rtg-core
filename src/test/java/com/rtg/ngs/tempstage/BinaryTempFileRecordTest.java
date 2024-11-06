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

import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class BinaryTempFileRecordTest extends TestCase {

  private BinaryTempFileRecord record() {
    final BinaryTempFileRecord bar = new BinaryTempFileRecord(false, true, false, false);
    bar.setAlignmentScore(1);
    bar.setCigarString("30=".getBytes());
    bar.setMdString("10A".getBytes());
    bar.setNumberMismatches(3);
    bar.setReadId(5);
    bar.setReferenceId(10);
    bar.setSamFlags((byte) 0);
    bar.setStartPosition(1000);
    return bar;
  }

  public void testBasics() throws IOException {
    final BinaryTempFileRecord bar = record();
    check(bar);
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (TempRecordWriterNio writer = new TempRecordWriterNio(baos)) {
      writer.writeRecord(bar);
      bar.setSentinelRecord();
      writer.writeRecord(bar);
    }
    assertTrue(bar.isSentinelRecord());
    try (TempRecordReaderNio dis = new TempRecordReaderNio(new ByteArrayInputStream(baos.toByteArray()), new TempRecordReader.RecordFactory(false, true, false, false))) {
      final BinaryTempFileRecord rec;
      assertNotNull(rec = dis.readRecord());
      check(rec);
      assertNull(dis.readRecord());
    }
  }

  private void check(BinaryTempFileRecord bar) {
    assertEquals(1, bar.getAlignmentScore());
    assertEquals("30=", new String(bar.getCigarString()));
    assertEquals("10A", new String(bar.getMdString()));
    assertEquals(3, bar.getNumberMismatches());
    assertEquals(5, bar.getReadId());
    assertEquals(10, bar.getReferenceId());
    assertEquals(0, bar.getSamFlags());
    assertEquals(1000, bar.getStartPosition());
    assertFalse(bar.isReadPaired());
    assertFalse(bar.isReverseStrand());
    assertFalse(bar.isSentinelRecord());
  }

  public void testIOExceptionBug1592() {
    final BinaryTempFileRecord bar = record();
    check(bar);
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final BufferedOutputStream bos = new BufferedOutputStream(baos) {
      @Override
      public synchronized void write(int b) throws IOException {
        throw new IOException("simulated out of disk space");
      }
      @Override
      public synchronized void write(byte[] b, int off, int len) throws IOException {
        throw new IOException("simulated out of disk space");
      }
    };
    final TempRecordWriterNio writer = new TempRecordWriterNio(bos);
    try {
      writer.writeRecord(bar);
      fail();
    } catch (final IOException e) {
      assertEquals("simulated out of disk space", e.getMessage());
    } finally {
      try {
        writer.close();
      } catch (final IOException e) {
        assertEquals("simulated out of disk space", e.getMessage());
      }
    }
  }

}
