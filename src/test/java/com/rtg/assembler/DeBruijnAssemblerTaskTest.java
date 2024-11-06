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
import java.util.Arrays;

import com.rtg.AbstractTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class DeBruijnAssemblerTaskTest extends AbstractTest {
  public void test() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File sequence = ReaderTestUtils.getDNADir(""
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCAATGAATCGAAAAAA" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCAATGAATCGAAAAAAGTGTG" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCAATGAATCGAAAAAAGTGTG" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCCATGAATCGAAAAAA" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCCATGAATCGAAAAAA" + StringUtils.LS
      );
      try {
        final MemoryPrintStream mps = new MemoryPrintStream();
        final MemoryPrintStream diag = new MemoryPrintStream();
        final DeBruijnParams params = DeBruijnParams.builder().inputFiles(Arrays.asList(sequence)).directory(tmpDir).kmerSize(6).create();
        final File fileList = new File(tmpDir, "fileList");
        FileUtils.stringToFile(sequence.toString() + StringUtils.LS, fileList);
        Diagnostic.setLogStream(diag.printStream());
        final DeBruijnAssemblerTask task = new DeBruijnAssemblerTask(params, mps.outputStream());
        task.exec();
//          System.err.println(FileUtils.fileToString(new File(tmpDir, "contigs")));
        TestUtils.containsAll(diag.toString()
          , "Maximum bubble length: 22"
          , "Tip threshold: 12"
        );
        final File contigs = new File(tmpDir, "contigs");
        assertTrue(contigs.exists());
        final File collapsed = new File(tmpDir, "collapsed");
        assertTrue(collapsed.exists());
        final File popped = new File(tmpDir, "popped");
        assertTrue(popped.exists());
      } finally {
        FileHelper.deleteAll(sequence);
      }
    }
  }
}
