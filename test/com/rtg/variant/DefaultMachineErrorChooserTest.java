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
