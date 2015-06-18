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
    UnmatedAugmenterCli.performAugmentation(mExistsSam, null, null, true, false, "fooaug", rgOut.outputStream());
    final String exp = FileHelper.resourceToString("com/rtg/sam/resources/mated.sam");
    final String res = FileUtils.fileToString(mExistsSam);

    final String expNoPg = exp.replaceAll("@PG.*\n", "");
    final String outStrNoPg = res.replaceAll("@PG.*\n", "");
    assertEquals(expNoPg, outStrNoPg);
    assertFalse(res.contains("@PG\tID:rtg"));
    TestUtils.containsAll(rgOut.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "rg1\t"
          );
  }

  public void testAugmentAsUnmated() throws IOException {
    FileHelper.resourceToFile("com/rtg/sam/resources/unmated.sam", mExistsSam);
    UnmatedAugmenterCli.performAugmentation(null, mExistsSam, null, true, false, "fooaug", null);
    final String exp = FileHelper.resourceToString("com/rtg/sam/resources/augmented.sam");
    final String res = FileUtils.fileToString(mExistsSam);

    final String expNoPg = TestUtils.stripSAMHeader(exp);
    final String outStrNoPg = TestUtils.stripSAMHeader(res);
    assertEquals(expNoPg, outStrNoPg);
  }

  public void testReadGroupStats() throws IOException {
    final MemoryPrintStream rgOut = new MemoryPrintStream();
    FileUtils.stringToFile(ReadGroupStatsCalculatorTest.RG_SAM, mExistsSam);
    UnmatedAugmenterCli.performAugmentation(mExistsSam, null, null, true, false, "fooaug", rgOut.outputStream());
    //flagrantly copied from ReadGroupStatsCalculatorTest
    TestUtils.containsAll(rgOut.toString(),
          "#CL",
          "#Version",
          ReadGroupStats.countsHeader(),
          "boo_pe\t"
          );
  }

  public void testAugmentFail() throws IOException {
    FileHelper.resourceToFile("com/rtg/sam/resources/unmated.sam", mExistsSam);
    assertTrue(mExistsSam.setReadOnly());
    try {
      try {
        UnmatedAugmenterCli.performAugmentation(mExistsSam, null, null, true, false, ".fooaug", null);
      } catch (final NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("rename"));
        assertTrue(e.getMessage().contains(mExistsSam.getName()));
        assertTrue(e.getMessage().contains("mated.fooaug.sam"));
      }
    } finally {
      mExistsSam.setWritable(true);
    }
  }

  public void testFileFilter() throws IOException {

    final File existsBamGz = new File(mDir, "mated.bam");
    assertTrue(existsBamGz.createNewFile());

    final File[] files = mDir.listFiles(new UnmatedAugmenterCli.SamBamGzMatedUnmatedFileFilter());
    assertEquals(3, files.length);
    TestUtils.containsAll(Arrays.toString(files), "unmated.sam.gz", "mated.sam", "mated.bam");
  }

  public final void testInitFlags() {
    checkHelp("Prepares SAM",
        "DIR");
  }

}
