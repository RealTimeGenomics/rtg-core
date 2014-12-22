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
package com.rtg.tabix;

import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class IndexUtilsTest extends TestCase {

  private static final String SAM_HEADER_CLIP = ""
    + "@HD" + TAB + "VN:1.4" + TAB + "SO:coordinate\n"
    + "@SQ" + TAB + "SN:t" + TAB + "LN:84\n"
    ;
  private static final String SAM_CLIP = SAM_HEADER_CLIP
    + "0" + TAB + "0" + TAB + "t" + TAB + "1" + TAB + "255" + TAB + "10S43M2S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAATCGCTAGGTTCGACTTGGTTAACAACAACGCCTGGGGCTTTTTGG" + TAB + "*\n"
    + "1" + TAB + "0" + TAB + "t" + TAB + "42" + TAB + "255" + TAB + "10S40M2S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAATTATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACGG" + TAB + "*\n";


  public void testCompress() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      ArrayList<File> files = new ArrayList<>();
      files.add(FileUtils.stringToFile(SAM_CLIP, new File(tmpDir, "sam.sam")));
      List<File> bzFiles = IndexUtils.ensureBlockCompressed(files);

      final File samFile = bzFiles.get(0);
      new TabixIndexer(samFile).saveSamIndex();
      assertTrue(samFile.exists());
      assertTrue(new File(samFile.toString() + ".tbi").exists());
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = new ExtractCli().mainInit(new String[] {samFile.getPath(), "--header"}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);

      assertEquals(SAM_CLIP, mps.toString());
    }
  }
}
