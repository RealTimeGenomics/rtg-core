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
package com.rtg.ml;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.util.Resources;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class ZeroRBuilderTest extends AbstractBuildClassifierTest {

  @Override
  BuildClassifier makeClassifier() {
    return new ZeroRBuilder();
  }

  public void testSaveLoad() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final ZeroRBuilder.ZeroRClassifier zero = new ZeroRBuilder.ZeroRClassifier(50, 79);
    zero.save(new DataOutputStream(mps.outputStream()), null);
    final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(mps.toByteArray()));
    assertEquals(MlPredictLoader.MlPredictType.ZERO_R.ordinal(), dis.readInt());
    final ZeroRBuilder.ZeroRClassifier zerob = new ZeroRBuilder.ZeroRClassifier(dis);
    assertEquals(zero, zerob);
  }

  private static final int MIN_VERSION = 1;
  public void testLoadVersionX() throws IOException {
    for (int i = MIN_VERSION; i <= ZeroRBuilder.ZeroRClassifier.SERIAL_VERSION; ++i) {
      checkLoadVersion(i);
    }
  }

  private void checkLoadVersion(int version) throws IOException {
    final InputStream is = Resources.getResourceAsStream("com/rtg/ml/resources/testZeroRVersion_" + version);
    try (final DataInputStream dis = new DataInputStream(is)) {
      assertEquals(MlPredictLoader.MlPredictType.ZERO_R.ordinal(), dis.readInt());
      ZeroRBuilder.ZeroRClassifier bs = new ZeroRBuilder.ZeroRClassifier(dis);
      TestUtils.containsAll(bs.toString(), "0R: 0.802");

      assertEquals(MlPredictLoader.MlPredictType.ZERO_R.ordinal(), dis.readInt());
      bs = new ZeroRBuilder.ZeroRClassifier(dis);
      TestUtils.containsAll(bs.toString(), "0R: 0.191");
    }
  }

  private static ZeroRBuilder.ZeroRClassifier createTestZeroR() {
    return new ZeroRBuilder.ZeroRClassifier(808, 199);
  }
  private static ZeroRBuilder.ZeroRClassifier createTestZeroR2() {
    return new ZeroRBuilder.ZeroRClassifier(19, 80);
  }

  public void testBadArgs() {
    try {
      new ZeroRBuilder.ZeroRClassifier(0, 0);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("0R undefined on no examples", e.getMessage());
    }
  }

  /**
   * Creates test serialized form for current serial version
   * @param args directory in which test data should be saved
   * @throws IOException it happens
   */
  public static void main(String[] args) throws IOException {
    final File dir = new File(args[0]);
    final File output = new File(dir, "testZeroRVersion_" + ZeroRBuilder.ZeroRClassifier.SERIAL_VERSION);
    try (DataOutputStream dos = new DataOutputStream(FileUtils.createOutputStream(output))) {
      createTestZeroR().save(dos, null);
      createTestZeroR2().save(dos, null);
    }
  }

}
