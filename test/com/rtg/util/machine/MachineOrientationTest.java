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

package com.rtg.util.machine;

import static com.rtg.sam.SharedSamConstants.OUT_SAM;
import static com.rtg.sam.SharedSamConstants.SAMHEADER1;
import static com.rtg.sam.SharedSamConstants.SAMHEADER2;
import static com.rtg.util.StringUtils.FS;
import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.sam.SamUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import junit.framework.TestCase;

/**
 */
public class MachineOrientationTest extends TestCase {


  public void test() {
    TestUtils.testEnum(MachineOrientation.class, "[FR, RF, TANDEM, ANY]");
  }
  private static final String TB = "\t";

  private static final String SAM_REC1 = ""
    + SAMHEADER1
    + "0" + TB +  "163" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "=" + TB + "11" + TB + "10" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC2 = ""
    + SAMHEADER1
   + "0" + TB + "83" + TB + "g1" + TB + "11" + TB + "255" + TB + "4M" + TB + "=" + TB + "3" + TB + "-10" + TB + "TTCA" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC3 = "" //
    + SAMHEADER2
    + "0" + TB +  "163" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "g2" + TB + "11" + TB + "10" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC4 = "" // reverse reverse second
    + SAMHEADER2
    + "0" + TB +  "179" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "=" + TB + "11" + TB + "10" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC5 = "" // Single end
    + SAMHEADER2
    + "0" + TB +  "16" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC6 = "" // mate unmapped
    + SAMHEADER1
    + "0" + TB +  "137" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "=" + TB + "11" + TB + "10" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC7 = "" // reverse reverse and different sequence
    + SAMHEADER2
    + "0" + TB +  "179" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "g2" + TB + "11" + TB + "10" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  private static final String SAM_REC8 = "" // reverse reverse
    + SAMHEADER2
    + "0" + TB +  "179" + TB + "g1" + TB +  "3" + TB + "255" + TB + "4M" + TB + "=" + TB + "11" + TB + "10" + TB + "ATCG" + TB + "````" + TB + "AS:i:1" + TB + "IH:i:1" + LS;

  public void testUnjumbleableIllumina() throws IOException {
    final MachineType illumina = MachineType.ILLUMINA_PE;
    final MachineOrientation orientationIllumina = illumina.orientation();
    SAMRecord rec = readOneSamRecord(SAM_REC1);
    assertTrue(orientationIllumina.firstOnTemplate(rec));
    assertTrue(orientationIllumina.hasValidMate(rec));

    SAMRecord rec2 = readOneSamRecord(SAM_REC2);
    assertFalse(orientationIllumina.firstOnTemplate(rec2));
    assertTrue(orientationIllumina.hasValidMate(rec2));

    SAMRecord rec3 = readOneSamRecord(SAM_REC3);
    assertFalse(orientationIllumina.hasValidMate(rec3));

    SAMRecord rec4 = readOneSamRecord(SAM_REC4);
    assertFalse(orientationIllumina.hasValidMate(rec4));

    SAMRecord rec5 = readOneSamRecord(SAM_REC5);
    assertFalse(orientationIllumina.hasValidMate(rec5));

    SAMRecord rec6 = readOneSamRecord(SAM_REC6);
    assertFalse(orientationIllumina.hasValidMate(rec6));

  }

  public void testUnjumbleableCg() throws IOException {
    final MachineType cg = MachineType.COMPLETE_GENOMICS;
    final MachineOrientation orientationCg = cg.orientation();
    SAMRecord rec = readOneSamRecord(SAM_REC8);
    assertTrue(orientationCg.hasValidMate(rec));

    SAMRecord rec7 = readOneSamRecord(SAM_REC7);
    assertFalse(orientationCg.hasValidMate(rec7));

    SAMRecord rec4 = readOneSamRecord(SAM_REC4);
    assertTrue(orientationCg.hasValidMate(rec4));
    assertTrue(orientationCg.firstOnTemplate(rec4));

    SAMRecord rec5 = readOneSamRecord(SAM_REC5);
    assertFalse(orientationCg.hasValidMate(rec5));

    SAMRecord rec6 = readOneSamRecord(SAM_REC6);
    assertFalse(orientationCg.hasValidMate(rec6));

  }

  private SAMRecord readOneSamRecord(String samstr) throws IOException {
    final File input = FileUtils.createTempDir("testcheck", "sv_in");
    try {
      FileUtils.stringToFile(samstr, new File(input, OUT_SAM));

      final String inn = input.getPath();

      final SAMRecord rec;
      File samfile = new File(inn + FS + OUT_SAM);
      try (SamReader reader = SamUtils.makeSamReader(FileUtils.createInputStream(samfile, false))) {
        final CloseableIterator<SAMRecord> iterator = reader.iterator();
        if (iterator.hasNext()) {
          rec = iterator.next();
        } else {
          rec = null;
        }

      }
      return rec;
    } finally {
      FileHelper.deleteAll(input);
    }
  }


  public void testJustReverse() {
    assertTrue(MachineOrientation.FR.orientationOkay(false, true, true, false));
    assertTrue(MachineOrientation.RF.orientationOkay(true, true, false, false));
    assertTrue(MachineOrientation.FR.orientationOkay(false, false, true, true));
    assertTrue(MachineOrientation.RF.orientationOkay(true, false, false, true));
    assertTrue(MachineOrientation.TANDEM.orientationOkay(true, false, true, true));
    assertTrue(MachineOrientation.TANDEM.orientationOkay(false, true, false, false));
    for (final boolean left : new boolean[] {false, true}) {
      for (final boolean right : new boolean[] {false, true}) {
        assertTrue(MachineOrientation.ANY.orientationOkay(left, true, right, false));
      }
    }
    assertFalse(MachineOrientation.FR.orientationOkay(true, true, false, false));
    assertFalse(MachineOrientation.FR.orientationOkay(true, false, true, true));
    assertFalse(MachineOrientation.FR.orientationOkay(false, true, false, false));
    assertFalse(MachineOrientation.RF.orientationOkay(false, false, true, true));
    assertFalse(MachineOrientation.RF.orientationOkay(true, true, true, false));
    assertFalse(MachineOrientation.RF.orientationOkay(false, true, false, false));
    assertFalse(MachineOrientation.TANDEM.orientationOkay(true, true, false, false));
    assertFalse(MachineOrientation.TANDEM.orientationOkay(false, false, true, true));
  }

  public void testPosAndReverse() {
    assertTrue(MachineOrientation.FR.orientationOkay(5, false, 8, true));
    assertTrue(MachineOrientation.RF.orientationOkay(10, true, 10, false));
    assertTrue(MachineOrientation.TANDEM.orientationOkay(10, true, 9, true));
    assertTrue(MachineOrientation.TANDEM.orientationOkay(9, false, 10, false));
    for (final boolean left : new boolean[] {false, true}) {
      for (final boolean right : new boolean[] {false, true}) {
        assertTrue(MachineOrientation.ANY.orientationOkay(0, left, 1, right));
      }
    }
    assertTrue(MachineOrientation.FR.orientationOkay(10, true, 5, false));
    assertFalse(MachineOrientation.FR.orientationOkay(4, true, 5, false));
    assertFalse(MachineOrientation.FR.orientationOkay(9, true, 10, true));
    assertFalse(MachineOrientation.FR.orientationOkay(5, false, 6, false));
    assertTrue(MachineOrientation.RF.orientationOkay(10, false, 5, true));
    assertFalse(MachineOrientation.RF.orientationOkay(10, false, 18, true));
    assertFalse(MachineOrientation.RF.orientationOkay(4, true, 3, true));
    assertFalse(MachineOrientation.RF.orientationOkay(2, false, 1, false));
    assertFalse(MachineOrientation.TANDEM.orientationOkay(1, true, 2, false));
    assertFalse(MachineOrientation.TANDEM.orientationOkay(4, false, 3, true));
    assertFalse(MachineOrientation.TANDEM.orientationOkay(9, true, 10, true));
    assertFalse(MachineOrientation.TANDEM.orientationOkay(10, false, 9, false));
  }
}
