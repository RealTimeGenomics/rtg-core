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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class DeBruijnAssemblerCliTest extends AbstractParamsCliTest<DeBruijnParams> {

  public void testFlags() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File sequence = ReaderTestUtils.getDNADir(">a" + StringUtils.LS + "ACGTACGTACGTACGTACGTACGTACGTAC" + StringUtils.LS);
      try {
        final File fileList = new File(tmpDir, "fileList");
        FileUtils.stringToFile(sequence.toString() + StringUtils.LS, fileList);

        checkHandleFlagsErr();
        final File output = new File(tmpDir, "bar");
        String err = checkHandleFlagsErr("-o", output.toString(), sequence.toString());
        assertTrue(err, err.contains(" You must provide a value for -k INT"));
        err = checkHandleFlagsErr("-o", output.toString(), "-k", "30");
        assertTrue(err, err.contains("No input files specified"));
        err = checkHandleFlagsErr(sequence.toString(), "-k", "30");
        assertTrue(err, err.contains("You must provide a value for -o DIR"));

        err = checkHandleFlagsErr("-o", output.toString(), "-I", fileList.toString(), "-k", "-30");
        assertTrue(err, err.contains("--kmer-size should be positive"));
        checkHandleFlags("-o", output.toString(), sequence.toString(), "-k", "30");

        checkMainInitOk("-o", output.toString(), "-I", fileList.toString(), "-k", "30");
        final CFlags flags =  new CFlags("foo", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
        DeBruijnAssemblerCli.initLocalFlags(flags);
        flags.setFlags("-o", output.toString(), "-I", fileList.toString(), "-k", "30", "-c", "42");
        final DeBruijnParams params = DeBruijnAssemblerCli.makeParamsLocal(flags);
        assertEquals(30, params.kmerSize());
        assertEquals(sequence, params.inputFiles().get(0));
        assertEquals(42, params.minHashFrequency());

      } finally {
        FileHelper.deleteAll(sequence);

      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testInitParams() {
    checkHelp("debruijn [OPTION]... -k INT -o DIR FILE+",
        "-k INT -o DIR -I FILE",
        "-k, --kmer-size=INT ", "kmer length to build graph nodes from",
        "-c, --minimum-kmer-frequency", "set minimum kmer frequency to retain"
        );
  }
  @Override
  protected ParamsCli<DeBruijnParams> getParamsCli() {
    return new DeBruijnAssemblerCli();
  }
}
