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
package com.rtg.tabix;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.regex.Pattern;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.Resources;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.gzip.GzipUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class IndexerCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new IndexerCli();
  }

  public void testUsage() {
    checkHelp("-f", "--format", "sam", "bam", "coveragetsv", "bed", "-I", "format of input", "input-list-file", "list of block compressed", "containing data", "TAB-delimited");
  }

  public void testFlags() throws IOException {
    final File f = FileUtils.createTempDir("foo", "rah");
    try {
      checkHandleFlagsOut("-f", "sam", f.getPath());
      checkHandleFlagsOut("-f", "bam", f.getPath());
      checkHandleFlagsOut("-f", "vcf", f.getPath());
      checkHandleFlagsOut("-f", "coveragetsv", f.getPath());
      checkHandleFlagsErr("-f", "hobbit", f.getPath());
      checkHandleFlagsErr(f.getPath());
    } finally {
      assertTrue(f.delete());
    }
  }

  public void testOperation() throws IOException {
    final File dir = FileUtils.createTempDir("indexercli", "test");
    try {
      final File file1 = FileHelper.resourceToFile("com/rtg/sam/resources/test.sam.gz", new File(dir, "test1.sam.gz"));
      final File file2 = FileHelper.resourceToFile("com/rtg/sam/resources/test.sam.gz", new File(dir, "test2.sam.gz"));
      final File file3 = FileHelper.resourceToFile("com/rtg/sam/resources/test.sam.gz", new File(dir, "test3.sam.gz"));
      final File file4;
      try (InputStream is = GzipUtils.createGzipInputStream(Resources.getResourceAsStream("com/rtg/sam/resources/test.sam.gz"))) {
        file4 = FileHelper.streamToFile(is, new File(dir, "test4.sam"));
      }
      try (MemoryPrintStream out = new MemoryPrintStream()) {
        final MemoryPrintStream err = new MemoryPrintStream();
        try {
          final int code = getCli().mainInit(new String[]{"-f", "sam", file1.getPath(), file2.getPath(), file3.getPath(), file4.getPath()}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 1, code);
          TestUtils.containsAll(StringUtils.grep(out.toString(), Pattern.quote(file1.getPath())), "Creating index for", "test1.sam.gz.tbi");
          TestUtils.containsAll(StringUtils.grep(out.toString(), Pattern.quote(file2.getPath())), "Creating index for", "test2.sam.gz.tbi");
          TestUtils.containsAll(StringUtils.grep(out.toString(), Pattern.quote(file3.getPath())), "Creating index for", "test3.sam.gz.tbi");
          assertEquals("Cannot create index for " + file4.getPath() + " as it is not in bgzip format." + StringUtils.LS, err.toString());
        } finally {
          err.close();
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testBamOperation() throws IOException {
    final File dir = FileUtils.createTempDir("indexercli", "test");
    try {
      final File file1 = FileHelper.resourceToFile("com/rtg/sam/resources/bam.bam", new File(dir, "bam1.bam"));
      final File file2 = FileHelper.resourceToFile("com/rtg/sam/resources/bam.bam", new File(dir, "bam2.bam"));
      final File file3 = FileHelper.resourceToFile("com/rtg/sam/resources/bam.bam", new File(dir, "bam3.bam"));
      final File file4 = FileHelper.resourceToFile("com/rtg/sam/resources/unmated.sam", new File(dir, "test4.sam"));
      try (MemoryPrintStream out = new MemoryPrintStream()) {
        final MemoryPrintStream err = new MemoryPrintStream();
        try {
          final int code = getCli().mainInit(new String[]{"-f", "bam", file1.getPath(), file2.getPath(), file3.getPath(), file4.getPath()}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 1, code);
          TestUtils.containsAll(StringUtils.grep(out.toString(), Pattern.quote(file1.getPath())), "Creating index for", "bam1.bam.bai");
          TestUtils.containsAll(StringUtils.grep(out.toString(), Pattern.quote(file2.getPath())), "Creating index for", "bam2.bam.bai");
          TestUtils.containsAll(StringUtils.grep(out.toString(), Pattern.quote(file3.getPath())), "Creating index for", "bam3.bam.bai");
          assertEquals("Cannot create index for " + file4.getPath() + " as it is not in bgzip format." + StringUtils.LS, err.toString());
        } finally {
          err.close();
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
