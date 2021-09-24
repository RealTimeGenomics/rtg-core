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

package com.rtg.util.store;

import java.io.File;
import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.io.TestDirectory;

/**
 */
public class StoreDirProxyTest extends AbstractStoreDirTest {

  @Override
  protected StoreDirectory getDirectory(File fileDir) throws IOException {
    final StoreDirProxy dir = new StoreDirProxy(fileDir);
    dir.integrity();
    return dir;
  }

  public void testBad() throws IOException {
    try (final TestDirectory tmp =  new TestDirectory()) {
      final File bad = new File(tmp, "bad");
      try {
        new StoreDirProxy(bad);
        fail();
      } catch (final IOException e) {
        final String message = e.getMessage();
        assertTrue(message.startsWith("Not a directory "));
        assertTrue(message.endsWith(StringUtils.FS + "bad"));
      }
    }
  }
}
