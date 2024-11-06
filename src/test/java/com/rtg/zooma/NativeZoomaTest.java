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
package com.rtg.zooma;

import java.io.File;
import java.io.IOException;

import com.rtg.AbstractTest;
import com.rtg.util.SpawnJvm;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.RandomDna;

/**
 * Test the corresponding class.
 */
public class NativeZoomaTest extends AbstractTest {

  /**
   * Allow capture of low level standard output in zooma.
   */
  public static final class TestNativeZooma {

    private TestNativeZooma() {
    }

    public static void main(final String[] args) throws IOException {
      try (final TestDirectory dir = new TestDirectory()) {
        final String reference = RandomDna.random(1000);
        final File referenceFile = new File(dir, "reference.fasta");
        final File indexFile = new File(dir, "reference.bin");
        FileUtils.stringToFile(">reference\n" + reference + "\n", referenceFile);
        final NativeZooma zooma = new NativeZooma();
        zooma.buildIndex(indexFile.getPath(), referenceFile.getPath(), null, null, 21, 1);
      }
    }
  }

  public void test() throws IOException {
    if (NativeZooma.isEnabled()) {
      final Process p = SpawnJvm.spawn(TestNativeZooma.class.getName());
      final String out = FileUtils.streamToString(p.getErrorStream());
      TestUtils.containsAll(out, "max nucleotides to hash: 1000", "saveindex time");
      assertEquals("", FileUtils.streamToString(p.getInputStream()));
      p.getOutputStream().close();
      p.getInputStream().close();
      p.getErrorStream().close();
    } else {
      System.err.println("skipping zooma test (no native lib)");
    }
  }
}
