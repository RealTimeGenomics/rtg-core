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
package com.rtg.ml;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class MlDataTypeTest extends TestCase {
  public void test() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    try (final DataOutputStream dos = new DataOutputStream(mps.outputStream())) {
      MlDataType.DOUBLE.save(5.5, dos);
      MlDataType.BOOLEAN.save(true, dos);
      MlDataType.INTEGER.save(556, dos);
      MlDataType.STRING.save("foorahgah", dos);
    }
    try (final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(mps.toByteArray()))) {
      assertEquals(5.5, MlDataType.DOUBLE.load(dis));
      assertEquals(true, MlDataType.BOOLEAN.load(dis));
      assertEquals(556, MlDataType.INTEGER.load(dis));
      assertEquals("foorahgah", MlDataType.STRING.load(dis));
    }
  }
}
