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
