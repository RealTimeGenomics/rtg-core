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
package com.rtg.ngs;


import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 * Test {@link CgMapCli}
 *
 */
public class CgMapFlagsValidatorTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CgMapCli();
  }

  public void testFailingMask2() throws Exception {
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(os);
    Diagnostic.setLogStream(ps);
    final File mainOut = FileUtils.createTempDir("cgmap", "test");
    try {
      final File left = new File(mainOut, "left");
      ReaderTestUtils.getReaderDNAFastqCG("@test\nacgt\n+\n####", left, PrereadArm.LEFT);
      final File right = new File(mainOut, "right");
      ReaderTestUtils.getReaderDNAFastqCG("@test\nacgt\n+\n####", right, PrereadArm.RIGHT);
      final File template = new File(mainOut, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgt", template, null);
      final File out = new File(mainOut, "out");
      final String[] args = {"-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath(), "--mask", "foobar"};
      final String err = checkHandleFlagsErr(args);
      assertTrue(err, err.contains("Invalid value \"foobar\" for flag --mask."));
    } finally {
      assertTrue(FileHelper.deleteAll(mainOut));
      Diagnostic.setLogStream();
    }
  }

  private static final String USAGE = ""
    + "Usage: rtg cgmap [OPTION]... -i SDF|FILE --mask STRING -o DIR -t SDF" + StringUtils.LS
    + StringUtils.LS
    + "Try '--help' for more information" + StringUtils.LS;

  public void testValidateFlags() throws Exception {
    final File mainOut = FileUtils.createTempDir("cgmap", "test");
    try {
      final File cgleft = new File(mainOut, "left");
      ReaderTestUtils.getReaderDNAFastqCG("@s1\ncgtccccgtacgtacgtatagatgctagcacgtac\n+\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n", cgleft, PrereadArm.LEFT);
      final File cgright = new File(mainOut, "right");
      ReaderTestUtils.getReaderDNAFastqCG("@s1\ngtacgtgctagcatcgtacgtacggacgtggtacg\n+\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n", cgright, PrereadArm.RIGHT);
      final File template = new File(mainOut, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgtacgtacgtacgtacgtacgtacgtgtc", template, null);
      final File testfile = new File(mainOut, "leftFile");
      FileUtils.stringToFile("b", testfile);
      final File out = new File(mainOut, "out");

      final String[] args = {"-i", "input",
        "--mask", "cg1",
        "-t", template.getPath(),
        "-o", out.getPath()
      };

      assertTrue(checkHandleFlagsErr(args).contains(USAGE));
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr(args), "The specified SDF, \"input\", does not exist");

      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-i", testfile.toString(), "-t", template.toString(), "-o", out.toString(), "--mask", "cg1"), "The specified file, \"" + testfile.toString() + "\", is not an SDF");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-i", mainOut.toString(), "-t", testfile.toString(), "-o", out.toString(), "--mask", "cg1"), "The specified file, \"" + testfile.toString() + "\", is not an SDF");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-i", cgleft.toString(), "-t", template.toString(), "-o", out.toString(), "--mask", "cg1"), "Inputfile not in paired end format.");

      checkHandleFlagsOut("-i", mainOut.toString(), "-t", template.toString(), "-o", out.toString(), "--mask", "cg1");
    } finally {
      assertTrue(FileHelper.deleteAll(mainOut));
    }
  }

}
