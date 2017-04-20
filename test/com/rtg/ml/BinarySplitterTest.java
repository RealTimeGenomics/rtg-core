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
      //assertEquals(version, bs.mCurrentVersion);
      TestUtils.containsAll(bs.toString(ds), "num6ers", String.valueOf(SPLIT_VALUE_1));

      bs = new BinarySplitter(dis, ds);
      //assertEquals(version, bs.mCurrentVersion);
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
