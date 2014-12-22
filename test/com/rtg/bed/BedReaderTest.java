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

package com.rtg.bed;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class BedReaderTest extends TestCase {

  private static final String BED_CONTENTS = ""
      + LS
      + "# COMMENT" + LS
      + "track adsfasdf" + LS
      + "browser adsfasdf" + LS
      + LS
      + "chr1" + TAB + "2" + TAB + "80" + TAB + "annotation" + LS
      + "chr2" + TAB + "3" + TAB + "40" + TAB + "annotation1" + TAB + "annotation2" + LS
      + LS
      + "chr3" + TAB + "7" + TAB + "90" + TAB + "annotation" + TAB + "annotation3" + LS;

  public void testBedReader() throws IOException {
    final String[] lines = BED_CONTENTS.split(LS);
    try (TestDirectory tempDir = new TestDirectory()) {
      final File bedFile = FileUtils.stringToFile(BED_CONTENTS, new File(tempDir, "test.bed"));
      try (BedReader reader = new BedReader(new BufferedReader(new FileReader(bedFile)))) {
        assertEquals(5, reader.getHeader().getHeaderLines().length);
        for (int i = 0; i < 5; i++) {
          assertEquals(lines[i], reader.getHeader().getHeaderLines()[i]);
        }
        assertTrue(reader.hasNext());
        assertEquals(lines[5], reader.next().toString());
        assertTrue(reader.hasNext());
        assertEquals(lines[6], reader.next().toString());
        assertTrue(reader.hasNext());
        assertEquals(lines[8], reader.next().toString());
        assertFalse(reader.hasNext());
      }
    }
  }
}
