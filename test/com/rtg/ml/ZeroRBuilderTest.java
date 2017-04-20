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
      final int type = dis.readInt();
      assertEquals(MlPredictLoader.MlPredictType.ZERO_R.ordinal(), type);
      final ZeroRBuilder.ZeroRClassifier bs = new ZeroRBuilder.ZeroRClassifier(dis);
      final StringBuilder str = bs.toString(new StringBuilder(), "", null);
      final String s = str.toString();
      assertTrue(s.contains("808/1007"));
    }
  }

  private static ZeroRBuilder.ZeroRClassifier createTestZeroR() {
    return new ZeroRBuilder.ZeroRClassifier(808, 199);
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
    }
  }

}
