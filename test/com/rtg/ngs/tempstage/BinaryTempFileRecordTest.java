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
