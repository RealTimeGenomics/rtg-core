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

import com.rtg.util.integrity.Exam;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class StoreFileResourceTest extends AbstractStoreTest {

  @Override
  protected StoreDirectory getDirectory(File tmpDir) {
    final StoreDirectory dir = new StoreDirResource("com/rtg/assembler/graph/io/resources");
    Exam.integrity(dir);
    return dir;
  }

  @Override
  public void testFile() throws IOException {
    final File unusedFile = new File("storeFileResourceDir");
    final StoreDirectory dir = getDirectory(unusedFile);
    assertFalse(unusedFile.exists());
    final StoreFile sf = dir.child("header1");
    assertEquals("header1", sf.name());
    try {
      sf.outputStream();
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }

    final String exp = FileHelper.resourceToString("com/rtg/assembler/graph/io/resources/header1");
    assertEquals(exp, sf.content());

    final InputStream is = sf.inputStream();
    assertEquals(exp, FileUtils.readerToString(new InputStreamReader(is)));
  }

}
