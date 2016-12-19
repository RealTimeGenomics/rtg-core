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
package com.rtg.variant.cnv;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Test class
 */
public class CnvProductCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CnvProductCli();
  }

  private File prepareTemplate(final File dir) throws IOException {
    return ReaderTestUtils.getDNADir(">t\na", dir);
  }

  public void test() {
    String res = checkHandleFlagsErr();
    final String exp = getCFlags().getUsageHeader();
    assertTrue(res.contains(exp));
    assertTrue(res.contains("You must provide a value for"));
    res = checkHandleFlagsErr("-i", new File("__googg").getPath());
    assertTrue(res.contains(exp));
    assertTrue(res.contains("You must provide a value for"));
  }

  public void test3() throws IOException {
    try (final TestDirectory dir = new TestDirectory("cnvcli")) {
      final File outputdir = new File(dir, "output");
      assertFalse(outputdir.exists());
      final String[] args = {
          "-i", new File(dir, "base").getPath(),
          "-j", new File(dir, "target").getPath(),
          "-t", prepareTemplate(dir).getPath(),
          "-o", outputdir.getPath(),
          "-b", "0",
      };
      final String x = checkHandleFlagsErr(args);
      assertTrue(x, x.contains("The bucket-size flag should be positive."));
    }
  }

  public void testMakeParams() throws Exception {
    try (final TestDirectory dir = new TestDirectory("cnvcli")) {
      final File[] baseExp = {
          new File(dir, "base"),
          new File(dir, "base2"),
          new File(dir, "base3")
      };
      for (final File f : baseExp) {
        assertTrue(f.createNewFile());
      }
      final File[] targetExp = {
          new File(dir, "target"),
          new File(dir, "target2"),
          new File(dir, "target3")
      };
      for (final File f : targetExp) {
        assertTrue(f.createNewFile());
      }
      final String[] args = {
          "-i", baseExp[0].getPath(),
          "-i", baseExp[1].getPath(),
          "-i", baseExp[2].getPath(),
          "-j", targetExp[0].getPath(),
          "-j", targetExp[1].getPath(),
          "-j", targetExp[2].getPath(),
          "-o", new File(dir, "output").getPath(),
          "-t", prepareTemplate(new File(dir, "template")).getPath(),
          "-b", "10",
          "-c", "1",
          "-m", "2",
          "-u", "3",
          "--Xallow-duplicate-start"
      };
      assertEquals("", checkHandleFlagsOut(args));
      final CnvProductParams params = ((CnvProductCli) mCli).makeParams();
      assertEquals(10, params.bucketSize());
      assertEquals(1, params.filterParams().maxAlignmentCount());
      assertEquals(new IntegerOrPercentage(2), params.filterParams().maxMatedAlignmentScore());
      assertEquals(new IntegerOrPercentage(3), params.filterParams().maxUnmatedAlignmentScore());
      assertEquals(false, params.filterStartPositions());
      assertEquals(new File(dir, "output"), params.directory());
      assertTrue(params.mappedBase().containsAll(Arrays.asList(baseExp)));
      assertTrue(params.mappedTarget().containsAll(Arrays.asList(targetExp)));
      assertEquals(3, params.mappedBase().size());
      assertEquals(3, params.mappedTarget().size());
    }
  }

  public void testInane() throws Exception {
    assertEquals("rtg", mCli.applicationName());
    assertEquals("cnv", mCli.moduleName());

    try (final TestDirectory dir = new TestDirectory("cnvcli")) {
      assertTrue(new File(dir, "base").createNewFile());
      assertTrue(new File(dir, "target").createNewFile());
      final String[] args = {
          "-i", new File(dir, "base").getPath(),
          "-t", prepareTemplate(new File(dir, "template")).getPath(),
          "-j", new File(dir, "target").getPath(),
          "-o", new File(dir, "output").getPath(),
      };
      assertEquals("", checkHandleFlagsOut(args));
      assertEquals(100, getCFlags().getValue(CnvProductCli.BUCKET_SIZE_FLAG));
      final CnvProductCli cli = (CnvProductCli) mCli;
      assertTrue(cli.task(cli.makeParams(), TestUtils.getNullOutputStream()) instanceof CnvProductTask);
      assertEquals(new File(dir, "output"), cli.outputDirectory());
    }
  }

}
