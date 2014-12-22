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

package com.rtg.visualization;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class SamHelperTest extends TestCase {

  public void testCorrectComplement() {
    assertEquals("ACGT", SamHelper.getCorrectComplement(true, "ACGT"));

    assertEquals("TTTT", SamHelper.getCorrectComplement(false, "AAAA"));
  }

  public void testValidAlignment() {
    final SAMRecord r = new SAMRecord(new SAMFileHeader());

    assertTrue(SamHelper.validAlignmentScore(r, 5, 4));
    r.setAttribute("AS", 3);

    assertTrue(SamHelper.validAlignmentScore(r, 5, 4));
    assertFalse(SamHelper.validAlignmentScore(r, 2, 4));

    final SAMRecord r2 = new SAMRecord(new SAMFileHeader());
    r2.setAttribute("AS", 3);
    r2.setMateUnmappedFlag(true);
    r2.setReadPairedFlag(true);

    assertTrue(SamHelper.validAlignmentScore(r2, 2, 4));
    assertFalse(SamHelper.validAlignmentScore(r2, 2, 2));
  }

  public void testValidIH() {
    final SAMRecord r = new SAMRecord(new SAMFileHeader());

    assertTrue(SamHelper.validIhScore(r, 5));

    r.setAttribute("IH", 3);

    assertTrue(SamHelper.validIhScore(r, 5));
    assertFalse(SamHelper.validIhScore(r, 2));
  }

  public void testloadAlignments() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("samhelper", "loadalignments");
    try {

      AviewTest.prepareData(f, AviewTest.SAM);

      final File[] alignments = {new File(f, AviewTest.ALIGNMENTS_BAM_FILE_NAME)};
      AviewParams p = new AviewParamsBuilder().alignments(alignments).region("gblah:5-6").create();
      try {
        SamHelper.loadAlignments(p);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("Unable to apply region \"gblah\" as the specified sequence name does not exist"));
      }
      p = new AviewParamsBuilder().alignments(alignments).region("g1:4-5").create();
      final ArrayList<SAMRecord> rec = SamHelper.loadAlignments(p);
      final int[] location = {3, 3, 5};
      assertEquals(location.length, rec.size());
      for (int i = 0; i < rec.size(); i++) {
          assertEquals(rec.get(i).getAlignmentStart(), location[i]);
      }
      alignments[0] = new File(f, AviewTest.ALIGNMENTS_SAM_FILE_NAME);
      p = new AviewParamsBuilder().alignments(alignments).region("g1:4-5").create();
      try {
        SamHelper.loadAlignments(p);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertEquals("File " + alignments[0].getPath() + " is not indexed", e.getMessage());
      }
    } finally {
      AviewTest.deleteBrokenBam(f);
    }

  }

  public void testSortAlignments() {

    final SAMRecord r = new SAMRecord(new SAMFileHeader());
    r.setAlignmentStart(5);
    r.setCigarString("2=");

    final SAMRecord r2 = new SAMRecord(new SAMFileHeader());
    r2.setAlignmentStart(1);
    r2.setCigarString("6=");

    final SAMRecord r3 = new SAMRecord(new SAMFileHeader());
    r3.setAlignmentStart(1);
    r3.setCigarString("15=");


    final ArrayList<SAMRecord> rec = new ArrayList<>();
    rec.add(r);
    rec.add(r2);
    rec.add(r3);

    SamHelper.sortAlignments(rec);

    assertEquals(r2, rec.get(0));
    assertEquals(r3, rec.get(1));
    assertEquals(r, rec.get(2));
  }

}
