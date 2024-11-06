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
