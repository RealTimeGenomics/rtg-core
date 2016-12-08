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
