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

package com.rtg.visualization;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class SamHelperTest extends TestCase {

  public void testCorrectComplement() {
    assertEquals("ACGT", SamHelper.getCorrectComplement(true, "ACGT"));

    assertEquals("TTTT", SamHelper.getCorrectComplement(false, "AAAA"));
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
        SamHelper.loadAlignments(p, null);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("Unable to apply region \"gblah\" as the specified sequence name does not exist"));
      }
      p = new AviewParamsBuilder().alignments(alignments).region("g1:4-5").create();
      final ArrayList<SAMRecord> rec = SamHelper.loadAlignments(p, null);
      final int[] location = {3, 3, 5};
      assertEquals(location.length, rec.size());
      for (int i = 0; i < rec.size(); ++i) {
          assertEquals(rec.get(i).getAlignmentStart(), location[i]);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(f));
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

    rec.sort(SamHelper.POSITION_COMP);

    assertEquals(r2, rec.get(0));
    assertEquals(r3, rec.get(1));
    assertEquals(r, rec.get(2));
  }

}
