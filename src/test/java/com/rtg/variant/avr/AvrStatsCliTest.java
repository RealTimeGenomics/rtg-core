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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class AvrStatsCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new AvrStatsCli();
  }

  public void test() throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final File model = FileHelper.resourceToFile("com/rtg/variant/avr/resources/mlmodel.avr", new File(dir, "mlmodel.avr"));
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = new AvrStatsCli().mainInit(new String[] {model.getPath()}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);
      mNano.check("avrstats.txt", StringUtils.grepMinusV(mps.toString(), "Location"));
    }
  }

}
