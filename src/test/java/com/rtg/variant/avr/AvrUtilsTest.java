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
package com.rtg.variant.avr;

import java.io.File;

import com.rtg.AbstractTest;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.TestDirectory;

public class AvrUtilsTest extends AbstractTest {

  public void testNoSuchAvrModel() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    flags.registerOptional(AvrUtils.AVR_MODEL_FILE_FLAG, File.class, "File", "AVR");
    flags.setFlags("--avr-model", "no-such-model-avr");
    try {
      AvrUtils.getAvrModel(flags, false);
      fail();
    } catch (final InvalidParamsException e) {
      // expected
    }
  }

  public void testAvrModelDirFailures() throws Exception {
    try (TestDirectory tmpDir = new TestDirectory()) {
      final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
      final String oldAvrModels = System.getProperty(AvrUtils.ENVIRONMENT_MODELS_DIR);
      try {
        final File f = new File(tmpDir, "blah");
        try {
          System.setProperty(AvrUtils.ENVIRONMENT_MODELS_DIR, f.getPath());
          AvrUtils.getAvrModel(flags, false);
          fail();
        } catch (final InvalidParamsException e) {
          assertFalse(f.exists());
          TestUtils.containsAllUnwrapped(e.getMessage(), "AVR models directory cannot be found or is not a directory: " + f.getPath());
        }
        assertTrue(f.createNewFile());
        try {
          System.setProperty(AvrUtils.ENVIRONMENT_MODELS_DIR, f.getPath());
          AvrUtils.getAvrModel(flags, false);
          fail();
        } catch (final InvalidParamsException e) {
          assertTrue(f.exists());
          TestUtils.containsAllUnwrapped(e.getMessage(), "AVR models directory cannot be found or is not a directory: " + f.getPath());
        }
      } finally {
        if (oldAvrModels == null) {
          System.clearProperty(AvrUtils.ENVIRONMENT_MODELS_DIR);
        } else {
          System.setProperty(AvrUtils.ENVIRONMENT_MODELS_DIR, oldAvrModels);
        }
      }
    }
  }
}
