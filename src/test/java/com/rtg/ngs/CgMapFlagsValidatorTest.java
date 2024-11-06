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
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-i", cgleft.toString(), "-t", template.toString(), "-o", out.toString(), "--mask", "cg1"), "Input file not in paired end format.");

      checkHandleFlagsOut("-i", mainOut.toString(), "-t", template.toString(), "-o", out.toString(), "--mask", "cg1");
    } finally {
      assertTrue(FileHelper.deleteAll(mainOut));
    }
  }

}
