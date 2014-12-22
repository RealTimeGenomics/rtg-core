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

import static com.rtg.util.StringUtils.FS;

import java.io.File;

import junit.framework.TestCase;

/**
 */
public class SdfFileUtilsTest extends TestCase {

  /**
   * Test method for {@link SdfFileUtils#sequenceDataFile(java.io.File, int)}.
   */
  public final void testSequenceDataFile() {
    final File dir = new File("dir");
    final File act = SdfFileUtils.sequenceDataFile(dir, 42);
    assertEquals("dir" + FS + "seqdata42", act.toString());
  }

  /**
   * Test method for {@link SdfFileUtils#sequencePointerFile(java.io.File, int)}.
   */
  public final void testSequencePointerFile() {
    final File dir = new File("dir");
    final File act = SdfFileUtils.sequencePointerFile(dir, 42);
    assertEquals("dir" + FS + "seqpointer42", act.toString());
  }

  /**
   * Test method for {@link SdfFileUtils#sequenceIndexFile(java.io.File)}.
   */
  public final void testSequenceIndexFile() {
    final File dir = new File("dir");
    final File act = SdfFileUtils.sequenceIndexFile(dir);
    assertEquals("dir" + FS + "sequenceIndex0", act.toString());
  }

  /**
   * Test method for {@link SdfFileUtils#labelDataFile(java.io.File, int)}.
   */
  public final void testLabelDataFile() {
    final File dir = new File("dir");
    final File act = SdfFileUtils.labelDataFile(dir, 42);
    assertEquals("dir" + FS + "namedata42", act.toString());
  }

  /**
   * Test method for {@link SdfFileUtils#labelPointerFile(java.io.File, int)}.
   */
  public final void testLabelPointerFile() {
    final File dir = new File("dir");
    final File act = SdfFileUtils.labelPointerFile(dir, 42);
    assertEquals("dir" + FS + "namepointer42", act.toString());
  }

  /**
   * Test method for {@link SdfFileUtils#labelIndexFile(java.io.File)}.
   */
  public final void testLabelIndexFile() {
    final File dir = new File("dir");
    final File act = SdfFileUtils.labelIndexFile(dir);
    assertEquals("dir" + FS + "nameIndex0", act.toString());
  }
}
