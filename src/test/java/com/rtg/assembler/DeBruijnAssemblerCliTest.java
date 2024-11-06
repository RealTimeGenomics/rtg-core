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
