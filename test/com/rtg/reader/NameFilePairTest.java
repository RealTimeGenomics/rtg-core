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

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class NameFilePairTest extends TestCase {

  public void testPairs() throws Exception {
    final File tmpDir = FileUtils.createTempDir("namefilepairtest", "blah");
    final File names = new File(tmpDir, "names");
    final File ptrs = new File(tmpDir, "ptrs");
    try {
      final NameFilePair nfp = new NameFilePair(names, ptrs, 13);

      assertTrue(nfp.canWriteName(1));
      assertFalse(nfp.canWriteName(13));

      nfp.writeName("woo!");

      assertTrue(names.exists());
      assertTrue(ptrs.exists());

      assertTrue(nfp.canWriteName(1));
      nfp.forceWriteName("feck");

      assertTrue(nfp.canWriteName(2));
      assertFalse(nfp.canWriteName(3));

      nfp.forceWriteName("BLLLAAA");

      assertFalse(nfp.canWriteName(1));

      nfp.writeName("fail");
      try {
        nfp.forceWriteName("");
        fail();
      } catch (final Exception e) {
        //expected
      } finally {
        nfp.close();
      }

      final String namesStr = FileUtils.fileToString(names);
      assertTrue(namesStr.contains("woo!"));
      assertTrue(namesStr.contains("feck"));
      assertTrue(namesStr.contains("BL"));
      assertFalse(namesStr.contains("BLL"));
      assertTrue(namesStr.contains("fail"));
      try (DataInputStream dis = new DataInputStream(new FileInputStream(ptrs))) {
        assertEquals(0, dis.readInt());
        assertEquals(5, dis.readInt());
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

}
