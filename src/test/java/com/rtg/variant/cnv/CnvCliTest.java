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
public class CnvCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CnvCli();
  }

  private File prepareTemplate(final File dir) throws IOException {
    return ReaderTestUtils.getDNADir(">t\na", dir);
  }

  public void test() {
    String res = checkHandleFlagsErr();
    assertTrue(res.contains("You must provide a value for"));
    res = checkHandleFlagsErr("-i", new File("__googg").getPath());
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
      final CnvProductParams params = ((CnvCli) mCli).makeParams();
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
      assertEquals(100, getCFlags().getValue(CnvCli.BUCKET_SIZE_FLAG));
      final CnvCli cli = (CnvCli) mCli;
      assertTrue(cli.task(cli.makeParams(), TestUtils.getNullOutputStream()) instanceof CnvProductTask);
      assertEquals(new File(dir, "output"), cli.outputDirectory());
    }
  }

}
