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

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class DefaultMachineErrorChooserTest extends TestCase {

  public void testDefault() {
    Diagnostic.setLogStream();
    final MachineErrorChooserInterface ch = new DefaultMachineErrorChooser();
    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    final VariantAlignmentRecord var = new VariantAlignmentRecord(r);
    assertNotNull(ch.machineErrors(var.getReadGroup(), var.isReadPaired()));
  }

  public void testSpecified() {
    Diagnostic.setLogStream();
    final MachineErrorParams me = MachineErrorParams.builder().create();
    final MachineErrorChooserInterface ch = new DefaultMachineErrorChooser(me);
    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    final VariantAlignmentRecord var = new VariantAlignmentRecord(r);
    assertTrue(me == ch.machineErrors(var.getReadGroup(), var.isReadPaired()));
  }

  public void testName() throws InvalidParamsException, IOException {
    Diagnostic.setLogStream();
    final MachineErrorChooserInterface ch = new DefaultMachineErrorChooser("complete");
    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    final VariantAlignmentRecord var = new VariantAlignmentRecord(r);
    assertNotNull(ch.machineErrors(var.getReadGroup(), var.isReadPaired()));
  }

  public void testNameBad() throws InvalidParamsException, IOException {
    Diagnostic.setLogStream();
    try {
      new DefaultMachineErrorChooser("xxx");
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertEquals("Invalid prior option \"xxx\"", e.getMessage());
    }
  }

  boolean mWarned = false;
  public void testPlatformWarnings() throws Exception {
    mWarned = false;
    final DiagnosticListener dl = new DiagnosticListener() {
      @Override
      public void handleDiagnosticEvent(DiagnosticEvent<?> event) {
        // warned only once
        assertFalse(mWarned);
        assertEquals(DefaultMachineErrorChooser.PLATFORM_WARNING, event.getMessage());
        mWarned = true;
      }

      @Override
      public void close() {
        // Nothing going on.
      }
    };
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    try {
      Diagnostic.addListener(dl);
      final MachineErrorChooserInterface chooser = new DefaultMachineErrorChooser("complete");
      final SAMFileHeader sfh = new SAMFileHeader();
      final SAMReadGroupRecord group1 = new  SAMReadGroupRecord("RG1");
      final SAMReadGroupRecord group2 = new  SAMReadGroupRecord("RG2");
      group1.setPlatform("A");
      group2.setPlatform("B");
      sfh.addReadGroup(group1);
      sfh.addReadGroup(group2);
      final SAMRecord r1 = new SAMRecord(sfh);
      r1.setAttribute("RG", group1.getId());
      final VariantAlignmentRecord var = new VariantAlignmentRecord(r1);
      chooser.machineErrors(var.getReadGroup(), var.isReadPaired());
      chooser.machineErrors(var.getReadGroup(), var.isReadPaired());
      assertFalse(mWarned);
      final SAMRecord r2 = new SAMRecord(sfh);
      r2.setAttribute("RG", group2.getId());
      final VariantAlignmentRecord var2 = new VariantAlignmentRecord(r2);
      chooser.machineErrors(var2.getReadGroup(), var2.isReadPaired());
      assertTrue(mWarned);
      chooser.machineErrors(var2.getReadGroup(), var2.isReadPaired());
    } finally {
      Diagnostic.removeListener(dl);
      Diagnostic.setLogStream();
    }
  }
}
