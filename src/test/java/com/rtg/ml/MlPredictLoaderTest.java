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
