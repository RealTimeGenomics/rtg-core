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
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.util.PortableRandom;
import com.rtg.util.Resources;
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
    TrainTestSplitTest.nukeData(data, 0.3);
    b.build(data);

    final PredictClassifier p = b.getClassifier();

    //System.err.println(p.toString(new StringBuilder(), ""));

    final SimpleEvaluation eval = new SimpleEvaluation();
    eval.evaluate(p, data);
    assertEquals(0.84, eval.accuracy(), 0.01);
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
      assertTrue(s.contains("12/57") && s.contains("99/112") && s.contains("808/809"));
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
