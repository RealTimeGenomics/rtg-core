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
public class BinarySplitterTest extends TestCase {

  private static final double SPLIT_VALUE_1 = 3213.6534;
  private static final double SPLIT_VALUE_2 = 3431.3322;

  public void testSaveLoad() throws IOException {
    final BinarySplitter bs = new BinarySplitter("random", 0, 44.44, MlDataType.DOUBLE);
    final MemoryPrintStream mps = new MemoryPrintStream();
    final DataOutputStream dos = new DataOutputStream(mps.outputStream());
    final Dataset data = new Dataset(new Attribute("random", MlDataType.DOUBLE));
    assertEquals("split: random <= 44.44", bs.toString(data));
    bs.save(dos, data);
    final BinarySplitter bs2 = new BinarySplitter(new DataInputStream(new ByteArrayInputStream(mps.toByteArray())), data);
    assertEquals(bs, bs2);
  }

  public void testSaveLoadMissing() throws IOException {
    final BinarySplitter bs = new BinarySplitter("random", 0, Double.NaN, MlDataType.DOUBLE);
    final MemoryPrintStream mps = new MemoryPrintStream();
    final DataOutputStream dos = new DataOutputStream(mps.outputStream());
    final Dataset data = new Dataset(new Attribute("random", MlDataType.DOUBLE));
    assertEquals("split: random == missing", bs.toString(data));
    bs.save(dos, data);
    final BinarySplitter bs2 = new BinarySplitter(new DataInputStream(new ByteArrayInputStream(mps.toByteArray())), data);
    assertEquals(bs, bs2);
  }

  private static final int MIN_VERSION = 1;
  public void testLoadVersionX() throws IOException {
    for (int i = MIN_VERSION; i <= BinarySplitter.SERIAL_VERSION; ++i) {
      checkLoadVersion(i);
    }
  }

  private void checkLoadVersion(int version) throws IOException {
    final Dataset ds = createTestDataset();
    final InputStream is = Resources.getResourceAsStream("com/rtg/ml/resources/testBinarySplitterVersion_" + version);
    try (final DataInputStream dis = new DataInputStream(is)) {
      BinarySplitter bs = new BinarySplitter(dis, ds);
      TestUtils.containsAll(bs.toString(ds), "num6ers", String.valueOf(SPLIT_VALUE_1));

      bs = new BinarySplitter(dis, ds);
      TestUtils.containsAll(bs.toString(ds), "num6ers", String.valueOf(SPLIT_VALUE_2));
    }
  }

  static Dataset createTestDataset() {
    return new Dataset(new Attribute("string with spaces and num6ers.", MlDataType.DOUBLE));
  }
  static BinarySplitter createTestSplitter() {
    return createTestSplitter(SPLIT_VALUE_1);
  }
  static BinarySplitter createTestSplitter(double splitValue) {
    final Dataset testDataset = createTestDataset();
    return new BinarySplitter(testDataset.getAttributes()[0].getName(), 0, splitValue, testDataset.getAttributes()[0].getDataType());
  }

  /**
   * Creates test serialized form for current serial version
   * @param args directory in which test data should be saved
   * @throws IOException it happens
   */
  public static void main(String[] args) throws IOException {
    final File dir = new File(args[0]);
    final File output = new File(dir, "testBinarySplitterVersion_" + BinarySplitter.SERIAL_VERSION);
    final Dataset ds = createTestDataset();
    try (DataOutputStream dos = new DataOutputStream(FileUtils.createOutputStream(output))) {
      createTestSplitter(SPLIT_VALUE_1).save(dos, ds);
      createTestSplitter(SPLIT_VALUE_2).save(dos, ds);
    }
  }
}
