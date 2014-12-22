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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;


/**
 */
public class BgZipTest extends AbstractCliTest {

  public void testHelp() {
    checkHelp("rtg bgzip",
        "Compress a file with block gzip",
        "--decompress", "file to (de)compress",
        "--stdout", "write on standard output, keep original files",
        "--force", "force overwrite of output file");
  }

  public void testCompress() throws Exception {
    final File dir = FileUtils.createTempDir("indexercli", "test");
    try {
      final File file1 = FileHelper.resourceToFile("com/rtg/sam/resources/test.sam", new File(dir, "test1.sam"));

      final BgZip cli = (BgZip) getCli();
      final MemoryPrintStream mps = new MemoryPrintStream();

      cli.mainInit(new String[] {file1.getPath()}, mps.outputStream(), mps.printStream());
      assertEquals(0, mps.toString().length());

      mps.reset();
      cli.mainInit(new String[] {"--decompress", "-c", new File(dir, "test1.sam.gz").getPath()}, mps.outputStream(), mps.printStream());

      assertTrue(mps.toString().length() > 1);
      assertEquals(FileHelper.resourceToString("com/rtg/sam/resources/test.sam"), mps.toString());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testCompressToStdOut() throws Exception {
    final File dir = FileUtils.createTempDir("indexercli", "test");
    try {
      final File file1 = FileHelper.resourceToFile("com/rtg/sam/resources/test.sam", new File(dir, "test1.sam"));

      final BgZip cli = (BgZip) getCli();
      final MemoryPrintStream mps = new MemoryPrintStream();

      cli.mainInit(new String[] {"--stdout", file1.getPath()}, mps.outputStream(), mps.printStream());
      assertTrue(mps.toString().length() > 1);
      assertTrue(file1.exists());

      final File bgzippedOutFile = new File(dir, "out.gz");
      FileUtils.byteArrayToFile(mps.outputStream().toByteArray(), bgzippedOutFile);

      mps.reset();
      cli.mainInit(new String[] {"-d", "-c", bgzippedOutFile.getPath()}, mps.outputStream(), mps.printStream());

      assertTrue(mps.toString().length() > 1);
      assertEquals(FileHelper.resourceToString("com/rtg/sam/resources/test.sam"), mps.toString());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testDecompress() throws Exception {
    final File dir = FileUtils.createTempDir("indexercli", "test");
    try {
      final File file1 = FileHelper.resourceToFile("com/rtg/sam/resources/test.sam.gz", new File(dir, "test1.sam.gz"));

      final BgZip cli = (BgZip) getCli();
      final MemoryPrintStream mps = new MemoryPrintStream();

      assertEquals(0, cli.mainInit(new String[] {"--decompress", file1.getPath()}, mps.outputStream(), mps.printStream()));

      assertEquals(FileHelper.resourceToString("com/rtg/sam/resources/test.sam"), FileUtils.fileToString(new File(dir, "test1.sam")));
      assertEquals(0, mps.toString().length());
      assertFalse(file1.exists());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testFailures() throws Exception {
    final File dir = FileHelper.createTempDirectory();
    try {
      final File file1 = new File(dir, "input");
      final BgZip cli = (BgZip) getCli();
      final MemoryPrintStream mps = new MemoryPrintStream();
      assertEquals(1, cli.mainInit(new String[] {"--decompress", file1.getPath()}, mps.outputStream(), mps.printStream()));

      assertEquals("Error: The specified file, \"" + file1.getPath() + "\" does not exist.", mps.toString().trim());

      FileUtils.stringToFile("", file1);
      final File fout = new File(dir, "input.gz");
      assertTrue(fout.createNewFile());

      mps.reset();
      assertEquals(1, cli.mainInit(new String[] {file1.getPath()}, mps.outputStream(), mps.printStream()));
      assertEquals("Error: Output file \"" + fout.getPath() + "\" already exists.", mps.toString().trim());
      assertEquals(0, cli.mainInit(new String[] {"-f", file1.getPath()}, mps.outputStream(), mps.printStream()));

      FileUtils.stringToFile("", file1);
      mps.reset();
      assertEquals(1, cli.mainInit(new String[] {"--decompress", file1.getPath()}, mps.outputStream(), mps.printStream()));
      assertEquals("Error: Input file not in GZIP format", mps.toString().trim());


    } finally {
      FileHelper.deleteAll(dir);
    }

  }

  @Override
  protected AbstractCli getCli() {
    return new BgZip();
  }
}
