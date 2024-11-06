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
import java.io.DataOutputStream;
import java.io.File;

import com.rtg.sam.SamBamConstants;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests corresponding class.
 */
public class SmartTempFileWriterTest extends TestCase {

  private File mDir;
  @Override
  public void setUp() throws Exception {
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    FileHelper.deleteAll(mDir);
    mDir = null;
  }

  public void testComparator() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final DataOutputStream dos = new DataOutputStream(mps.outputStream());
    final SmartTempFileWriter writer = new SmartTempFileWriter(dos, new PairedTempFileRecordComparator(), 8);
    assertEquals(0, writer.getMaxCapacityUsed());

    final BinaryTempFileRecord samrec1 = makeRecord(1, 1, true, 6, false);
    final BinaryTempFileRecord samrec2 = makeRecord(1, 1, false, 6, true);
    final BinaryTempFileRecord samrec3 = makeRecord(1, 1, false, 6, false);
    final BinaryTempFileRecord samrec4 = makeRecord(1, 1, false, 9, false);
    final BinaryTempFileRecord samrec5 = makeRecord(1, 2, false, 6, false);
    final BinaryTempFileRecord samrec6 = makeRecord(5, 1, false, 6, false);
    final BinaryTempFileRecord samrec7 = makeRecord(5, 1, false, 6, false);
    final BinaryTempFileRecord samrec8 = makeRecord(5, 1, false, 6, true);

    assertTrue(writer.addAlignmentHandleDuplicates(samrec4));
    assertTrue(writer.addAlignmentHandleDuplicates(samrec5));
    assertTrue(writer.addAlignmentHandleDuplicates(samrec2));
    assertTrue(writer.addAlignmentHandleDuplicates(samrec3));
    assertTrue(writer.addAlignmentHandleDuplicates(samrec8));
    assertTrue(writer.addAlignmentHandleDuplicates(samrec6));
    assertTrue(writer.addAlignmentHandleDuplicates(samrec1));
    assertFalse(writer.addAlignmentHandleDuplicates(samrec6));
    assertFalse(writer.addAlignmentHandleDuplicates(samrec6));
    assertFalse(writer.addAlignmentHandleDuplicates(samrec7));
    assertFalse(writer.addAlignmentHandleDuplicates(samrec7));
    writer.close();

    assertEquals(4, writer.getDuplicateCount());
    //check the file contents I guess..
    final TempRecordReader dis = new TempRecordReaderNio(new ByteArrayInputStream(mps.toByteArray()), new TempRecordReader.RecordFactory(true, false, false, false));

    BinaryTempFileRecord r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 1, (byte) 17, 1, 6);
    r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 1, (byte) 33, 1, 6);
    r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 1, (byte) 1, 1, 6);
    r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 1, (byte) 1, 1, 9);
    r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 2, (byte) 1, 1, 6);
    r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 1, (byte) 33, 5, 6);
    r = dis.readRecord();
    assertNotNull(r);
    checkRecord(r, 1, (byte) 1, 5, 6);

    assertNull(dis.readRecord());

  }

  public void testOverflow() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final DataOutputStream dos = new DataOutputStream(mps.outputStream());
    final SmartTempFileWriter writer = new SmartTempFileWriter(dos, new PairedTempFileRecordComparator(), 8);
    assertEquals(0, writer.getMaxCapacityUsed());
    for (int i = 0; i <= 10; ++i) {
      writer.addAlignmentHandleDuplicates(makeRecord(i + 1, 1, false, i + 5, false));
    }

    final BinaryTempFileRecord samrecoverflow = makeRecord(1, 1, false, 6, false);
    try {
      writer.addAlignmentHandleDuplicates(samrecoverflow);
      fail();
    } catch (final IllegalStateException e) {
      // Expected
    }

    assertEquals(10, writer.getMaxCapacityUsed());

  }

  void checkRecord(BinaryTempFileRecord r, int readId, byte flags, int startPosition, int mateAlignStart) {
    assertEquals(readId, r.getReadId());
    assertEquals(0, r.getReferenceId());
    assertEquals(flags, r.getSamFlags());
    assertEquals(startPosition, r.getStartPosition());
    assertEquals(mateAlignStart, r.getMatePosition());

  }

  BinaryTempFileRecord makeRecord(int alignStart, int readId, boolean negStrand, int mateAlignStart, boolean mateNegStrand) {
    final BinaryTempFileRecord rec = new BinaryTempFileRecord(true, false, false, false);
    rec.setStartPosition(alignStart);
    rec.setReferenceId(0);
    rec.setReadId(readId);
    rec.setMatePosition(mateAlignStart);
    rec.setCigarString(new byte[0]);
    rec.setMdString(new byte[0]);
    int samFlags = 0;
    if (negStrand) {
      samFlags |= SamBamConstants.SAM_READ_IS_REVERSE;
    }
      samFlags |= SamBamConstants.SAM_READ_IS_PAIRED;
    if (mateNegStrand) {
      samFlags |= SamBamConstants.SAM_MATE_IS_REVERSE;
    }
    rec.setSamFlags((byte) samFlags);
    return rec;
  }
}
