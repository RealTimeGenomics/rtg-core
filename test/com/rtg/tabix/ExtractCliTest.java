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

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;

/**
 * Test class
 */
public class ExtractCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ExtractCli();
  }

  public void testHelp() {
    checkHelp("rtg extract"
        , "Extract records from an indexed genomic position data file."
        , "FILE the indexed block compressed genome position data file to extract"
        , "STRING+ the range to display. The format is one of <template_name>, <template_name>:start-end or <template_name>:start+length. May be specified 0 or more times Reporting"
        , "--header", "print out header also"
        , "-h,", "--help", "print help on command-line flag usage"
        );
  }

  public void testSomeFlags() {
    checkHandleFlagsErr();
    checkHandleFlagsOut("foo");
    checkHandleFlagsOut("foo", "bar");
    checkHandleFlagsOut("foo", "bar", "bang");
  }

  public void testNormal() throws IOException {
    try (TestDirectory dir = new TestDirectory("extractcli")) {
      final File tabix = new File(dir, "snp_only.vcf.gz.tbi");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz.tbi", tabix);
      final File input = new File(dir, "snp_only.vcf.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz", input);
      try (MemoryPrintStream out = new MemoryPrintStream()) {
        try (MemoryPrintStream err = new MemoryPrintStream()) {
          int code = getCli().mainInit(new String[]{input.getPath(), "simulatedSequence19:500-1000"}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 0, code);
          mNano.check("extract500-1000", out.toString());

          out.reset();
          err.reset();
          code = getCli().mainInit(new String[]{input.getPath(), "simulatedSequence19:1000-5000"}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 0, code);
          mNano.check("extract1000-5000", out.toString());

          out.reset();
          err.reset();
          code = getCli().mainInit(new String[]{input.getPath(), "simulatedSequence19:500-5000"}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 0, code);
          mNano.check("extract500-5000", out.toString());

          out.reset();
          err.reset();
          code = getCli().mainInit(new String[]{input.getPath(), "simulatedSequence19:500-1000", "simulatedSequence19:1000-5000"}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 0, code);
          mNano.check("extract500-1000-5000", out.toString());

          out.reset();
          err.reset();
          code = getCli().mainInit(new String[]{input.getPath()}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 0, code);
          mNano.check("extract-norestrict", out.toString());

          out.reset();
          err.reset();
          code = getCli().mainInit(new String[]{input.getPath(), "--header-only"}, out.outputStream(), err.printStream());
          assertEquals(err.toString(), 0, code);
          mNano.check("extract-header-only", out.toString());
        }
      }
    }
  }

  public void testErrors() throws IOException {
    try (TestDirectory dir = new TestDirectory("extractcli")) {
      File input = FileUtils.stringToFile("stuff", new File(dir, "afile"));
      try (MemoryPrintStream out = new MemoryPrintStream()) {
        try (MemoryPrintStream err = new MemoryPrintStream()) {
          int code = getCli().mainInit(new String[]{input.getPath()}, out.outputStream(), err.printStream());
          assertEquals(1, code);
          assertTrue(err.toString().contains("bgzip format"));
          assertTrue(err.toString().contains(input.getPath()));
          input = BgzipFileHelper.bytesToBgzipFile("stuff".getBytes(), input);
          out.reset();
          err.reset();
          code = getCli().mainInit(new String[]{input.getPath()}, out.outputStream(), err.printStream());
          assertEquals(1, code);
          assertTrue(err.toString().contains("Index not found"));
          assertTrue(err.toString().contains(input.getPath()));
          assertTrue(err.toString().contains(input.getPath() + TabixIndexer.TABIX_EXTENSION));
        }
      }
    }
  }

}
