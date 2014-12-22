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
package com.rtg.util.io;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * JUnit tests for the corresponding class.
 *
 */
public class PartitionTest extends TestCase {

  public void testNormalBehaviour() throws Exception {
    final File main = FileUtils.createTempDir("partition", "test");
    try {
      final ArrayList<File> list = new ArrayList<>();
      final File[] f = new File[10];
      final StringBuilder sb = new StringBuilder();
      for (int k = 0; k < f.length; k++) {
        sb.append("X");
        f[k] = new File(main, String.valueOf(k));
        FileUtils.stringToFile(sb.toString(), f[k]);
        list.add(f[k]);
      }
      final List<List<File>> bins = Partition.partition(3, f);
      assertEquals(3, bins.size());
      assertEquals(f[0], bins.get(2).get(3));
      assertEquals(f[1], bins.get(2).get(2));
      assertEquals(f[2], bins.get(1).get(2));
      assertEquals(f[3], bins.get(0).get(2));
      assertEquals(f[4], bins.get(0).get(1));
      assertEquals(f[5], bins.get(1).get(1));
      assertEquals(f[6], bins.get(2).get(1));
      assertEquals(f[7], bins.get(2).get(0));
      assertEquals(f[8], bins.get(1).get(0));
      assertEquals(f[9], bins.get(0).get(0));

      final List<List<File>> bins2 = Partition.partition(3, list);
      assertEquals(3, bins2.size());
      assertEquals(f[0], bins2.get(2).get(3));
      assertEquals(f[1], bins2.get(2).get(2));
      assertEquals(f[2], bins2.get(1).get(2));
      assertEquals(f[3], bins2.get(0).get(2));
      assertEquals(f[4], bins2.get(0).get(1));
      assertEquals(f[5], bins2.get(1).get(1));
      assertEquals(f[6], bins2.get(2).get(1));
      assertEquals(f[7], bins2.get(2).get(0));
      assertEquals(f[8], bins2.get(1).get(0));
      assertEquals(f[9], bins2.get(0).get(0));
} finally {
      assertTrue(FileHelper.deleteAll(main));
    }
  }

}
