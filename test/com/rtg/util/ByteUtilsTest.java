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
package com.rtg.util;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 */
public class ByteUtilsTest extends TestCase {

  public void testConstants() throws IOException {
    final OutputStream out = new ByteArrayOutputStream();
    ByteUtils.writeLn(out);
    out.write(ByteUtils.FS_BYTE);
    out.write(ByteUtils.TAB_BYTE);
    assertEquals(StringUtils.LS + StringUtils.FS + "\t", out.toString());
  }
}
