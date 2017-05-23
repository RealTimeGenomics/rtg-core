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
