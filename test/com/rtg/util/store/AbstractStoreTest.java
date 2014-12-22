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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;

import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractStoreTest extends TestCase {

  protected abstract StoreDirectory getDirectory(File fileDir) throws IOException;

  public void testFile() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final String str = "foobar";
      final StoreDirectory dir = getDirectory(tmpDir);
      final StoreFile sf = dir.child("child1");
      assertEquals("child1", sf.name());
      final OutputStream os = sf.outputStream();
      os.write(str.getBytes());
      os.close();

      assertEquals(str, sf.content());

      final InputStream is = sf.inputStream();
      assertEquals(str, FileUtils.readerToString(new InputStreamReader(is)));
    }
  }
}
