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
import java.io.IOException;

import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SdfUtilsTest extends TestCase {

  public void test() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("sdfutils", "unittest");
    try {
      try {
        SdfUtils.validateHasNames(tempDir);
        fail();
      } catch (NoTalkbackSlimException e) {
        TestUtils.containsAll(e.getMessage(), "Unable to find file: " + new File(tempDir, "mainIndex") + " (", ") part of SDF: " + tempDir);
      }
      final File templ = ReaderTestUtils.getDNADir(">seq1\nactg", new File(tempDir, "template"), false, false, false);
      try {
        SdfUtils.validateHasNames(templ);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertEquals("SDF: " + templ + " has no name data", e.getMessage());
      }
      SequenceParams params = SequenceParams.builder().directory(templ).loadNames(true).mode(SequenceMode.UNIDIRECTIONAL).create();
      try {
        SdfUtils.validateNoDuplicates(params, false);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertEquals("SDF: " + templ + " has no name data", e.getMessage());
      }
      assertTrue(FileHelper.deleteAll(templ));
      ReaderTestUtils.getDNADir(">seq1\nactg\n>seq1\ngtca", templ);
      SdfUtils.validateHasNames(templ);
      params = SequenceParams.builder().directory(templ).loadNames(true).mode(SequenceMode.UNIDIRECTIONAL).create();
      try {
        SdfUtils.validateNoDuplicates(params, false);
        fail();
      } catch (NoTalkbackSlimException e) {
        assertEquals("Duplicate sequence names detected in SDF: " + templ, e.getMessage());
      } finally {
        params.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }
}
