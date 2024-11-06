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
package com.rtg.variant.sv;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class UnmatedAugmenterCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new UnmatedAugmenterCli();
  }

  private File mDir = null;
  private File mExistsSam = null;
  private File mExistsSamGz = null;
  private File mExistsNormal = null;
  private File mList = null;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileUtils.createTempDir("unmated", "augmenter");
    mExistsNormal = new File(mDir, "anormalfile");
    mExistsSam = new File(mDir, "mated.sam");
    mExistsSamGz = new File(mDir, "unmated.sam.gz");
    assertTrue(mExistsNormal.createNewFile());
    assertTrue(mExistsSam.createNewFile());
    assertTrue(mExistsSamGz.createNewFile());
    mList = new File(mDir, "filelist");
    FileUtils.stringToFile(mExistsNormal.getPath() + StringUtils.LS
            + mExistsSam.getPath() + StringUtils.LS
            + mExistsSamGz.getPath() + StringUtils.LS
            , mList);
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
    mExistsNormal = null;
    mExistsSam = null;
    mExistsSamGz = null;
    mList = null;
  }

  public void testFlagsBad() {
    checkHandleFlagsErr();
    checkHandleFlagsErr("-k");
    checkHandleFlagsErr("-k", "-s", "foobarsuffix");
    checkHandleFlagsErr(mDir.getPath(), "-s", "foobarsuffix");
    checkHandleFlagsErr(mExistsNormal.getPath());

  }

  public void testFlagsOk() {
    checkHandleFlagsOut(mDir.getPath());
    checkHandleFlagsOut(mDir.getPath(), "-k");
    checkHandleFlagsOut(mDir.getPath(), "-k", "-s", "foobarsuffix");
  }

  public void testGetOutputFile() throws IOException {
    File f = UnmatedAugmenterCli.getOutputFile(mExistsNormal, false, false, ".aug");
    assertEquals(mExistsNormal.getName() + ".aug", f.getName());
    f = UnmatedAugmenterCli.getOutputFile(mExistsSam, false, false, ".aug");
    assertEquals("mated.aug.sam", f.getName());
    f = UnmatedAugmenterCli.getOutputFile(mExistsSamGz, false, false, ".aug");
    assertEquals("unmated.aug.sam.gz", f.getName());
    assertTrue(f.createNewFile());
    try {
      f = UnmatedAugmenterCli.getOutputFile(mExistsSamGz, false, false, ".aug");
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertTrue(e.getMessage().contains("File: "));
      assertTrue(e.getMessage().contains("exists"));
      assertTrue(e.getMessage().contains(f.getName()));
    }
    final File oldF = f;
    f = UnmatedAugmenterCli.getOutputFile(mExistsSamGz, false, true, ".aug");
    assertEquals(oldF, f);
    f = UnmatedAugmenterCli.getOutputFile(mExistsSamGz, true, false, ".aug");
    assertFalse(f.equals(oldF));
    assertEquals(mDir, f.getParentFile());
    assertTrue(f.getName().contains(".temp"));
  }

  public void testAugmentMated() throws IOException {
    final MemoryPrintStream rgOut = new MemoryPrintStream();
    FileHelper.resourceToFile("com/rtg/sam/resources/mated.sam", mExistsSam);
    FileHelper.resourceToFile("com/rtg/sam/resources/unmated.sam", mExistsSamGz);
    UnmatedAugmenterCli.augmentSplitFiles(mExistsSam, mExistsSamGz, null, true, false, "fooaug", rgOut.outputStream());

    final String res = FileUtils.fileToString(mExistsSam);
    assertFalse(res.contains("@PG\tID:rtg"));
    final String expNoPg = FileHelper.resourceToString("com/rtg/sam/resources/mated.sam").replaceAll("@PG.*\n", "");
    final String outStrNoPg = res.replaceAll("@PG.*\n", "");
    assertEquals(expNoPg, outStrNoPg);
    TestUtils.containsAll(rgOut.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "rg1\t"
          );

    assertEquals(TestUtils.stripSAMHeader(FileHelper.resourceToString("com/rtg/sam/resources/augmented.sam")),
      TestUtils.stripSAMHeader(FileHelper.gzFileToString(mExistsSamGz)));
  }

  static final String EMPTY_SAM = ""
    + "@HD\tVN:1.0\tSO:coordinate" + StringUtils.LS
    + "@SQ\tSN:simulatedSequence1\tLN:400" + StringUtils.LS;

  public void testReadGroupStats() throws IOException {
    final MemoryPrintStream rgOut = new MemoryPrintStream();
    FileUtils.stringToFile(ReadGroupStatsCalculatorTest.RG_SAM, mExistsSam);
    FileHelper.stringToGzFile(EMPTY_SAM, mExistsSamGz);
    UnmatedAugmenterCli.augmentSplitFiles(mExistsSam, mExistsSamGz, null, true, false, "fooaug", rgOut.outputStream());
    //flagrantly copied from ReadGroupStatsCalculatorTest
    TestUtils.containsAll(rgOut.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "boo_pe\t"
          );
  }

  public void testAugmentFail() throws IOException {
    FileHelper.resourceToFile("com/rtg/sam/resources/mated.sam", mExistsSam);
    FileHelper.resourceToGzFile("com/rtg/sam/resources/unmated.sam", mExistsSamGz);
    assertTrue(mExistsSam.setReadOnly());
    try {
      try {
        UnmatedAugmenterCli.augmentSplitFiles(mExistsSam, mExistsSamGz, null, true, false, ".fooaug", new MemoryPrintStream().outputStream());
      } catch (final NoTalkbackSlimException e) {
        assertTrue(e.getMessage(), e.getMessage().contains("rename"));
        assertTrue(e.getMessage().contains(mExistsSam.getName()));
        assertTrue(e.getMessage().contains("mated.fooaug.sam"));
      }
    } finally {
      mExistsSam.setWritable(true);
    }
  }

  public void testFileFilters() throws IOException {
    assertTrue(new File(mDir, "mated.bam").createNewFile());
    assertTrue(new File(mDir, "unmapped.bam").createNewFile());

    File[] files = mDir.listFiles(new UnmatedAugmenterCli.MatedFileFilter());
    assertNotNull(files);
    assertEquals(2, files.length);
    TestUtils.containsAll(Arrays.toString(files), "mated.sam", "mated.bam");

    files = mDir.listFiles(new UnmatedAugmenterCli.UnmatedFileFilter());
    assertNotNull(files);
    assertEquals(1, files.length);
    TestUtils.containsAll(Arrays.toString(files), "unmated.sam.gz");

    files = mDir.listFiles(new UnmatedAugmenterCli.UnmappedFileFilter());
    assertNotNull(files);
    assertEquals(1, files.length);
    TestUtils.containsAll(Arrays.toString(files), "unmapped.bam");
  }

  public final void testInitFlags() {
    checkHelp("Prepares SAM",
        "DIR");
  }

}
