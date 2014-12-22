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
import java.io.DataOutputStream;
import java.io.File;

import com.rtg.sam.SamBamConstants;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

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
    for (int i = 0; i <= 10; i++) {
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


  public static Test suite() {
    return new TestSuite(SmartTempFileWriterTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
