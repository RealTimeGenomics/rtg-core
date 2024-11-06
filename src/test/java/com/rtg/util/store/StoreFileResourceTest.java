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
