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
package com.rtg.variant.sv;




import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.OutputParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;
import com.rtg.variant.sv.SvToolParams.SvToolParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class SvToolParamsTest extends TestCase {

  public void testOmnes() {
    new TestParams(SvToolParams.class, SvToolParamsBuilder.class).check();
  }

  public void testDefaults() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("svtoolparams", "test");
    try {
      final SvToolParams def = SvToolParams.builder().mapped(Collections.emptyList()).outputParams(new OutputParams(tempDir, false, false)).genome(new MockReaderParams(1, 1, SequenceMode.BIDIRECTIONAL)).create();
      assertEquals(def.binSize(), 0);
      assertEquals(def.stepSize(), 0);
      assertEquals(def.fineStepSize(), 0);
      assertFalse(def.heterozygous());
      assertFalse(def.outputSimple());
      assertNull(def.correctionsFile());
      try (OutputStream bayStream = def.bayesianStream()) {
        bayStream.write("blah blah".getBytes());
      }
      final File bayFile = new File(tempDir, "sv_bayesian.tsv");
      assertTrue(bayFile.isFile());
      assertEquals("blah blah", FileUtils.fileToString(bayFile));
      TestUtils.containsAll(def.toString()
          , "SvToolParams "
          , " output-simple=" + false
          , " heterozygous=" + false
          , " bin-size=" + 0
          , " step-size=" + 0
          , " fine-step-size=" + 0
          , " correctionsFile=" + null
          );
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testBuilder() {
    final SvToolParamsBuilder builder = SvToolParams.builder();
    assertEquals(builder, builder.binSize(1));
    assertEquals(builder, builder.stepSize(1));
    assertEquals(builder, builder.fineStepSize(1));
    assertEquals(builder, builder.heterozygous(true));
    assertEquals(builder, builder.outputSimple(true));
    assertEquals(builder, builder.correctionsFile(null));
  }
}
