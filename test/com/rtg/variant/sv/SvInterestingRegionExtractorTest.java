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
package com.rtg.variant.sv;

import java.io.File;

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class SvInterestingRegionExtractorTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testIRE() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("sviret")) {
      final File out = new File(tmpDir, "blah.bed.gz");
      final File in = new File(tmpDir, "in.txt");
      FileHelper.resourceToFile("com/rtg/variant/sv/resources/sv_int_region.txt", in);

      SvInterestingRegionExtractor.main(new String[] {"-i", in.getPath(), "-o", out.getPath()});

      final String s = FileHelper.gzFileToString(out);

      TestUtils.containsAll(s,
          "#chr\tstart\tend\tareas\tmaxscore\taverage",
          "chr1\t10\t4750\t3\t171.7091\t85.0291",
          "chr1\t7340\t9970\t2\t98.1045\t33.6914",
          "chr1\t39600\t0\t1\t186.4456\t159.4061",
          "chr2\t220493770\t220495380\t2\t52.3520\t22.2873"
          );

      final File tabixFile = new File(tmpDir, "blah.bed.gz.tbi");
      assertTrue(tabixFile.exists());
    }
  }

  public void testIREMerge() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("sviret")) {
      final File out = new File(tmpDir, "blah.bed");
      final File in = new File(tmpDir, "in.txt");
      FileHelper.resourceToFile("com/rtg/variant/sv/resources/sv_int_region.txt", in);

      SvInterestingRegionExtractor.main(new String[] {"-i", in.getPath(), "-o", out.getPath(), "-m", "3000"});

      final String s = FileUtils.fileToString(out);

      TestUtils.containsAll(s,
          "#chr\tstart\tend\tareas\tmaxscore\taverage",
          "chr1\t10\t9970\t6\t171.7091\t67.9166",
          "chr1\t39600\t0\t1\t186.4456\t159.4061",
          "chr2\t220493770\t220495380\t2\t52.3520\t22.2873"
      );
      final File tabixFile = new File(tmpDir, "blah.bed.tbi");
      assertFalse(tabixFile.exists());
      final File tabixFile2 = new File(tmpDir, "blah.bed.gz.tbi");
      assertFalse(tabixFile2.exists());
    }
  }
}
