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
package com.rtg.variant;

import java.io.ByteArrayInputStream;
import java.io.IOException;

import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 * Test class
 */
public class ReadGroupMachineErrorChooserTest extends TestCase {

  private static final String SAM = ""
    + "@RG\tSM:foo\tID:sounique\tPL:ILLUMINA" + StringUtils.LS
    + "@RG\tSM:foo\tID:uniquer\tPL:COMPLETE" + StringUtils.LS
    + "@RG\tSM:foo\tID:somereads\tPL:LS454" + StringUtils.LS
    + "@RG\tSM:foo\tID:paired\tPL:LS454\tPI:500" + StringUtils.LS
    + "@RG\tSM:foo\tID:differentyetthesame\tPL:ILLUMINA" + StringUtils.LS;

  public void test() throws Exception {
    final SAMFileHeader header = SamUtils.makeSamReader(new ByteArrayInputStream(SAM.getBytes())).getFileHeader();
    final ReadGroupMachineErrorChooser m = new ReadGroupMachineErrorChooser(header);
    final SAMRecord s = new SAMRecord(header);
    s.setAttribute("RG", "sounique");
    VariantAlignmentRecord var = new VariantAlignmentRecord(s);
    assertEquals(MachineType.ILLUMINA_PE, m.machineErrors(var.getReadGroup(), var.isReadPaired()).machineType());
    s.setAttribute("RG", "uniquer");
    var = new VariantAlignmentRecord(s);
    assertEquals(MachineType.COMPLETE_GENOMICS, m.machineErrors(var.getReadGroup(), var.isReadPaired()).machineType());
    s.setAttribute("RG", "somereads");
    var = new VariantAlignmentRecord(s);
    assertEquals(MachineType.FOURFIVEFOUR_SE, m.machineErrors(var.getReadGroup(), var.isReadPaired()).machineType());
    s.setAttribute("RG", "paired");
    var = new VariantAlignmentRecord(s);
    assertEquals(MachineType.FOURFIVEFOUR_PE, m.machineErrors(var.getReadGroup(), var.isReadPaired()).machineType());
    s.setAttribute("RG", "differentyetthesame");
    var = new VariantAlignmentRecord(s);
    assertEquals(MachineType.ILLUMINA_PE, m.machineErrors(var.getReadGroup(), var.isReadPaired()).machineType());
  }

  private static final String SAM_BAD = ""
    + "@RG\tSM:foo\tID:sounique" + StringUtils.LS;

  public void testMissingPlatform() throws IOException {
    Diagnostic.setLogStream();
    final SAMFileHeader header = SamUtils.makeSamReader(new ByteArrayInputStream(SAM_BAD.getBytes())).getFileHeader();
    try {
      new ReadGroupMachineErrorChooser(header);
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertTrue(e.getMessage().startsWith("Read group: sounique has no specified platform"));
    }
  }

  public void testSamRecord() throws Exception {
    Diagnostic.setLogStream();
    final SAMFileHeader header = SamUtils.makeSamReader(new ByteArrayInputStream(SAM.getBytes())).getFileHeader();
    final MachineErrorChooserInterface m = new ReadGroupMachineErrorChooser(header);
    final SAMReadGroupRecord rgr = new SAMReadGroupRecord("somethingelse");
    rgr.setPlatform("IONTORRENT");
    header.addReadGroup(rgr);
    final SAMRecord s = new SAMRecord(header);
    try {
      final VariantAlignmentRecord var = new VariantAlignmentRecord(s);
      m.machineErrors(var.getReadGroup(), var.isReadPaired());
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertEquals("Sam record had no read group attribute, but header read groups were supplied.", e.getMessage());
    }
    s.setAttribute("RG", "somethingelse");
    rgr.setPlatform("BLAH");
    try {
      final VariantAlignmentRecord var = new VariantAlignmentRecord(s);
      m.machineErrors(var.getReadGroup(), var.isReadPaired());
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertEquals("Sam record referenced read group \"somethingelse\" which was not found in the header.", e.getMessage());
    }
  }
}
