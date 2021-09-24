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

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class UnfilteredBinaryTempFileRecordTest extends TestCase {

  public void testIt() throws Exception {

    final BinaryTempFileRecord bupar = new BinaryTempFileRecord(true, false, false, true);
    bupar.setReferenceId(8);
    bupar.setCigarString(new byte[0]);
    bupar.setMdString(new byte[0]);
    assertFalse(bupar.isUnfilteredMated());

    final MemoryPrintStream mps = new MemoryPrintStream();
    final TempRecordWriterNio writer = new TempRecordWriterNio(mps.outputStream());
    writer.writeRecord(bupar);

    BinaryTempFileRecord bupar2;
    TempRecordReaderNio reader = new TempRecordReaderNio(new ByteArrayInputStream(mps.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, true));
    assertNotNull(bupar2 = reader.readRecord());

    assertEquals(8, bupar2.getReferenceId());
    assertFalse(bupar2.isUnfilteredMated());

    mps.reset();

    bupar.setUnfilteredMated(true);
    assertTrue(bupar.isUnfilteredMated());
    writer.writeRecord(bupar);
    reader = new TempRecordReaderNio(new ByteArrayInputStream(mps.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, true));
    assertNotNull(bupar2 = reader.readRecord());
    assertTrue(bupar2.isUnfilteredMated());

    mps.reset();

    bupar.setSentinelRecord();
    writer.writeRecord(bupar);
    reader = new TempRecordReaderNio(new ByteArrayInputStream(mps.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, true));
    assertNull(reader.readRecord());
  }

}
