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
