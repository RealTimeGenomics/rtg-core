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
package com.rtg.graph;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;

/**
 */
public class RocPlotCliTest extends AbstractCliTest {
  @Override
  protected AbstractCli getCli() {
    return new RocPlotCli();
  }

  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      checkHandleFlagsErr();
      final File f = new File(dir, "moo");
      checkHandleFlagsErr(f.getPath());
      assertTrue(f.createNewFile());
      checkHandleFlagsOut(f.getPath());
      checkHandleFlagsOut("--curve", f.getPath() + "=Monkey");
      final File png = new File(dir, "png.png");
      checkHandleFlagsOut(f.getPath(), "--png", png.getPath());
      assertTrue(png.createNewFile());
      checkHandleFlagsErr(f.getPath(), "--png", png.getPath());

      final String html = checkMainInitOk("--curve", f.getPath() + "=Monkey", "--Xhtml");
      TestUtils.containsAll(html,
              String.format("<param name=\"data1\" value=\"%s\">", f.getPath()),
              String.format("<param name=\"name1\" value=\"%s\">", "Monkey")
      );
      final File f2 = new File(dir, "oink");
      assertTrue(f2.createNewFile());
      checkHandleFlagsOut("--curve", f.getPath() + "=Monkey", f2.getPath());

      final String html2 = checkMainInitOk("--curve", f.getPath() + "=Monkey", f2.getPath(), "--Xapplet");
      TestUtils.containsAll(html2,
              String.format("<param name=\"data1\" value=\"%s\">", f.getPath()),
              String.format("<param name=\"name1\" value=\"%s\">", "Monkey"),
              String.format("<param name=\"data2\" value=\"%s\">", f2.getPath()),
              String.format("<param name=\"name2\" value=\"%s\">", f2.getParentFile().getName())
      );
    }
  }
}
