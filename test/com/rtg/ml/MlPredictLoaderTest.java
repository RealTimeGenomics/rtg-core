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
import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class MlPredictLoaderTest extends TestCase {
  public void testZeroRLoad() throws IOException {
    final ZeroRBuilder.ZeroRClassifier zero = new ZeroRBuilder.ZeroRClassifier(54, 93);
    checkSaveLoad(zero, null);
  }

  public void testBinaryTreeLoad() throws IOException {
    final Dataset data = new Dataset(new Attribute("booboo", MlDataType.INTEGER));
    final PredictClassifier bin = new BinaryTreeClassifier(new BinarySplitter(data.getAttributes()[0].getName(), 0, 5, data.getAttributes()[0].getDataType()), new ZeroRBuilder.ZeroRClassifier(33, 77), new ZeroRBuilder.ZeroRClassifier(99, 55), 0.44);
    checkSaveLoad(bin, data);
  }

  public void testBaggingLoad() throws IOException {
    final PredictClassifier bag = new BaggedClassifier(new ZeroRBuilder.ZeroRClassifier(11, 41), new ZeroRBuilder.ZeroRClassifier(999, 80), new ZeroRBuilder.ZeroRClassifier(1337, 6006));
    checkSaveLoad(bag, null);
  }

  static void checkSaveLoad(PredictClassifier bin, Dataset data) throws IOException {
    final StringBuilder exp = bin.toString(new StringBuilder(), "", data);
    final MemoryPrintStream mps = new MemoryPrintStream();
    bin.save(new DataOutputStream(mps.outputStream()), data);
    mps.outputStream().close();
    final PredictClassifier clazz = MlPredictLoader.loadPredictClassifier(new ByteArrayInputStream(mps.toByteArray()), data);
    final StringBuilder actual = clazz.toString(new StringBuilder(), "", data);
    assertEquals(exp.toString(), actual.toString());
  }
  public void testCorrupt() {
    try {
      checkCorruptLoad(1, (byte) -1);
      fail();
    } catch (IOException e) {
      TestUtils.containsAll(e.getMessage(), "Prediction type out of range");

    }
    try {
      checkCorruptLoad(BinarySplitter.SERIAL_VERSION == 1 ? 24 : 18, (byte) -1);
      fail();
    } catch (IOException e) {
      TestUtils.containsAll(e.getMessage(), "Learning attribute out of range");
    }
  }
  public void checkCorruptLoad(int pos, byte value) throws IOException {
    final Dataset data = new Dataset(new Attribute("booboo", MlDataType.INTEGER));
    final PredictClassifier bin = new BinaryTreeClassifier(new BinarySplitter(data.getAttributes()[0].getName(), 0, 5, data.getAttributes()[0].getDataType()), new ZeroRBuilder.ZeroRClassifier(33, 77), new ZeroRBuilder.ZeroRClassifier(99, 55), 0.44);
    final MemoryPrintStream mps = new MemoryPrintStream();
    bin.save(new DataOutputStream(mps.outputStream()), data);
    mps.outputStream().close();
    final byte[] array = mps.toByteArray();
    array[pos] = value;
    MlPredictLoader.loadPredictClassifier(new ByteArrayInputStream(array), data);
  }
}
