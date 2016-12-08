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

import com.rtg.util.Resources;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class BinaryTreeClassifierTest extends TestCase {

  public void testTreePredictorManual() {
    final BinaryTreeClassifier node2 = new BinaryTreeClassifier(new BinarySplitter("att1", 1, 1.0, MlDataType.BOOLEAN),
        new ZeroRBuilder.ZeroRClassifier(0.7),
        new ZeroRBuilder.ZeroRClassifier(1.0),
        0.2
    );
    final BinaryTreeClassifier classifier = new BinaryTreeClassifier(new BinarySplitter("att0", 0, Math.PI, MlDataType.DOUBLE),
        new ZeroRBuilder.ZeroRClassifier(0.9),
        node2,
        0.3
    );

    //System.err.println(classifier.toString(new StringBuilder(), ""));

    final double tolerance = 0.00001;
    //double, boolean
    assertEquals(0.9, classifier.predict(new double[] {3.0, 1.0}));
    assertEquals(0.9, classifier.predict(new double[] {Math.PI, 1.0}));
    assertEquals(0.9, classifier.predict(new double[] {3.0, 0.0}));
    assertEquals(0.9, classifier.predict(new double[] {3.0, Double.NaN}));

    assertEquals(0.7, classifier.predict(new double[] {3.5, 1.0}));
    assertEquals(1.0, classifier.predict(new double[] {3.5, 0.0}));
    assertEquals(0.94, classifier.predict(new double[] {3.5, Double.NaN}), tolerance);

    assertEquals(0.76, classifier.predict(new double[] {Double.NaN, 1.0}), tolerance);
    assertEquals(0.97, classifier.predict(new double[] {Double.NaN, 0.0}), tolerance);
    assertEquals(0.928, classifier.predict(new double[] {Double.NaN, Double.NaN}), tolerance);
  }


  public void testSaveLoad() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final Dataset data = new Dataset(new Attribute("lala", MlDataType.DOUBLE));
    final BinaryTreeClassifier bin = new BinaryTreeClassifier(new BinarySplitter(data.getAttributes()[0].getName(), 0, 33.3, data.getAttributes()[0].getDataType()), new ZeroRBuilder.ZeroRClassifier(50, 79), new ZeroRBuilder.ZeroRClassifier(90, 123), 0.90);
    bin.save(new DataOutputStream(mps.outputStream()), data);
    final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(mps.toByteArray()));
    assertEquals(MlPredictLoader.MlPredictType.BINARY_TREE.ordinal(), dis.readInt());
    final BinaryTreeClassifier binb = new BinaryTreeClassifier(dis, data);
    assertEquals(bin, binb);
  }

  private static final int MIN_VERSION = 1;
  public void testLoadVersionX() throws IOException {
    for (int i = MIN_VERSION; i <= BinaryTreeClassifier.SERIAL_VERSION; ++i) {
      checkLoadVersion(i);
    }
  }

  private void checkLoadVersion(int version) throws IOException {
    final InputStream is = Resources.getResourceAsStream("com/rtg/ml/resources/testBinaryTreeVersion_" + version);
    try (final DataInputStream dis = new DataInputStream(is)) {
      final int type = dis.readInt();
      assertEquals(MlPredictLoader.MlPredictType.BINARY_TREE.ordinal(), type);
      final BinaryTreeClassifier bs = new BinaryTreeClassifier(dis, BinarySplitterTest.createTestDataset());
      assertEquals(version, bs.mCurrentVersion);
      final StringBuilder str = bs.toString(new StringBuilder(), "", BinarySplitterTest.createTestDataset());
      final String s = str.toString();
      assertTrue(s.contains("54/97") && s.contains("90/190"));
    }
  }

  private static BinaryTreeClassifier createTestTree() {
    return new BinaryTreeClassifier(BinarySplitterTest.createTestSplitter(), new ZeroRBuilder.ZeroRClassifier(54, 43), new ZeroRBuilder.ZeroRClassifier(90, 100), 0.43);
  }

  /**
   * Creates test serialized form for current serial version
   * @param args directory in which test data should be saved
   * @throws IOException it happens
   */
  public static void main(String[] args) throws IOException {
    final File dir = new File(args[0]);
    final File output = new File(dir, "testBinaryTreeVersion_" + BinaryTreeClassifier.SERIAL_VERSION);
    try (DataOutputStream dos = new DataOutputStream(FileUtils.createOutputStream(output))) {
      createTestTree().save(dos, BinarySplitterTest.createTestDataset());
    }
  }
}
