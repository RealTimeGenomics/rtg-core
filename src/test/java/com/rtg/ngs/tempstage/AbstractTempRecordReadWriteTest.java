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
import java.io.InputStream;
import java.io.OutputStream;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractTempRecordReadWriteTest extends TestCase {

  protected abstract TempRecordWriter getWriter(OutputStream out, TempRecordReader.RecordFactory fact);
  protected abstract TempRecordReader getReader(InputStream in, TempRecordReader.RecordFactory fact);

  public void testReadWrite() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final TempRecordWriter wr = getWriter(baos, new TempRecordReader.RecordFactory(true, true, true, true));
    for (int i = 0; i < 30; ++i) {
      final BinaryTempFileRecord rec = new BinaryTempFileRecord(true, true, true, true);
      rec.setStartPosition(i);
      rec.setReadDeltaString("foo".getBytes());
      rec.setCgReadString("bar".getBytes());
      rec.setSuperCigarString("bang".getBytes());
      rec.setAlignmentScore(1);
      rec.setCigarString("baz".getBytes());
      rec.setMdString("moo".getBytes());
      rec.setNumberMismatches(1);
      rec.setReadId(7);
      rec.setReferenceId(10);
      rec.setSamFlags((byte) 30);
      rec.setComboScore(4);
      rec.setMatePosition(90);
      rec.setTemplateLength(100);
      wr.writeRecord(rec);
    }
    wr.close();
    final ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
    final TempRecordReader reader = getReader(bais, new TempRecordReader.RecordFactory(true, true, true, true));
    for (int i = 0; i < 30; ++i) {
      final BinaryTempFileRecord rec = reader.readRecord();
      assertEquals(rec.getStartPosition(), i);
      assertEquals("foo", new String(rec.getReadDeltaString()));
      assertEquals("bar", new String(rec.getCgReadString()));
      assertEquals("bang", new String(rec.getSuperCigarString()));
      assertEquals(1, rec.getAlignmentScore());
      assertEquals("baz", new String(rec.getCigarString()));
      assertEquals("moo", new String(rec.getMdString()));
      assertEquals(7, rec.getReadId());
      assertEquals(10, rec.getReferenceId());
      assertEquals(30, rec.getSamFlags());
      assertEquals(4, rec.getComboScore());
      assertEquals(90, rec.getMatePosition());
      assertEquals(100, rec.getTemplateLength());
    }
    reader.close();
  }
}
