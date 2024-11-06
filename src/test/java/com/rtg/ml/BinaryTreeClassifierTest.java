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

import junit.framework.TestCase;

/**
 */
public class BinaryTreeClassifierTest extends TestCase {

  public void testTreePredictorManual() {
    final BinaryTreeClassifier node2 = new BinaryTreeClassifier(new BinarySplitter("att1", 1, 1.0, MlDataType.BOOLEAN),
        new ZeroRBuilder.ZeroRClassifier(7, 3),
        new ZeroRBuilder.ZeroRClassifier(10, 0),
        0.2
    );
    final BinaryTreeClassifier classifier = new BinaryTreeClassifier(new BinarySplitter("att0", 0, Math.PI, MlDataType.DOUBLE),
        new ZeroRBuilder.ZeroRClassifier(9, 1),
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
      final StringBuilder str = bs.toString(new StringBuilder(), "", BinarySplitterTest.createTestDataset());
      final String s = str.toString();
      TestUtils.containsAll(s, "0R: 0.556", "0R: 0.473");
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
