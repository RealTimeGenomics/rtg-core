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

package com.rtg.reader;

import java.io.File;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SamToFastqTest extends AbstractNanoTest {

  public void test() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File inputFile = FileHelper.resourceToFile("com/rtg/reader/resources/input.sam", new File(dir, "input.sam"));
      final File outputPrefix = new File(dir, "output");
      final MainResult mr = MainResult.run(new SamToFastq(), "-i", inputFile.getPath(), "-o", outputPrefix.getPath(), "-Z");
      assertEquals(0, mr.rc());
      assertEquals("", mr.err());
      final File leftFastqOut = new File(dir, "output_1.fastq");
      final File rightFastqOut = new File(dir, "output_2.fastq");
      mNano.check("expected_1.fastq", FileUtils.fileToString(leftFastqOut));
      mNano.check("expected_2.fastq", FileUtils.fileToString(rightFastqOut));
    }
  }
}
