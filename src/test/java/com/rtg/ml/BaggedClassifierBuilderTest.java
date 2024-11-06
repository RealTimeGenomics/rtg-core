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

import com.rtg.util.PortableRandom;
import com.rtg.util.Resources;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class BaggedClassifierBuilderTest extends AbstractBuildClassifierTest {

  @Override
  BuildClassifier makeClassifier() {
    return new BaggedClassifierBuilder();
  }

  public void testOob() {
    final BuildClassifier b = makeClassifier();
    final Dataset data = TrainTestSplitTest.makeSimpleDataset(100, 100);
    b.build(data);
    final double oob = ((BaggedClassifierBuilder) b).getOutOfBagAccuracy();
    assertTrue(oob >= 0 && oob <= 1.0);
    //System.err.println(b.getClassifier().toString(new StringBuilder(), ""));
  }

  public void testCircleTree() {
    final BuildClassifier b = makeClassifier();

    final Dataset data = TrainTestSplitTest.makeCircleDataset(new PortableRandom(42), 100, 200);
    data.injectMissing(0.3);
    b.build(data);

    final PredictClassifier p = b.getClassifier();

    //System.err.println(p.toString(new StringBuilder(), ""));

    final SimpleEvaluation eval = new SimpleEvaluation();
    eval.evaluate(p, data);
    assertEquals(0.94, eval.accuracy(), 0.01);
  }

  public void testSaveLoadPredict() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final BaggedClassifier bag = new BaggedClassifier(new ZeroRBuilder.ZeroRClassifier(50, 79), new ZeroRBuilder.ZeroRClassifier(90, 123), new ZeroRBuilder.ZeroRClassifier(111, 222));
    bag.save(new DataOutputStream(mps.outputStream()), null);
    final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(mps.toByteArray()));
    assertEquals(MlPredictLoader.MlPredictType.BAGGED.ordinal(), dis.readInt());
    final BaggedClassifier bagb = new BaggedClassifier(dis, null);
    assertEquals(bag, bagb);
  }

  private static final int MIN_VERSION = 1;
  public void testLoadVersionX() throws IOException {
    for (int i = MIN_VERSION; i <= BaggedClassifier.SERIAL_VERSION; ++i) {
      checkLoadVersion(i);
    }
  }

  private void checkLoadVersion(int version) throws IOException {
    final InputStream is = Resources.getResourceAsStream("com/rtg/ml/resources/testBaggedVersion_" + version);
    try (final DataInputStream dis = new DataInputStream(is)) {
      final int type = dis.readInt();
      assertEquals(MlPredictLoader.MlPredictType.BAGGED.ordinal(), type);
      final BaggedClassifier bs = new BaggedClassifier(dis, null);
      final StringBuilder str = bs.toString(new StringBuilder(), "", null);
      final String s = str.toString();
      TestUtils.containsAll(s, "classifier [1]", "0R: 0.210", "0R: 0.883", "0R: 0.998");
    }
  }

  private static BaggedClassifier createTestBagged() {
    return new BaggedClassifier(new ZeroRBuilder.ZeroRClassifier(12, 45), new ZeroRBuilder.ZeroRClassifier(99, 13), new ZeroRBuilder.ZeroRClassifier(808, 1));
  }

  /**
   * Creates test serialized form for current serial version
   * @param args directory in which test data should be saved
   * @throws IOException it happens
   */
  public static void main(String[] args) throws IOException {
    final File dir = new File(args[0]);
    final File output = new File(dir, "testBaggedVersion_" + BaggedClassifier.SERIAL_VERSION);
    try (DataOutputStream dos = new DataOutputStream(FileUtils.createOutputStream(output))) {
      createTestBagged().save(dos, null);
    }
  }

}
