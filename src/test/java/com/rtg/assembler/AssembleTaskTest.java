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
import java.util.Collections;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class AssembleTaskTest extends TestCase {
  public void testAssembleTask() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final String[] s = new String[30];
      for (int i = 0; i < s.length; ++i) {
        s[i] = "ACGGGACATACGTGTGAGATACGATAGCACAGGACGTGATGACGTCCCGTCG";
      }
      final File input =  ReaderTestUtils.getDNADir(ReaderTestUtils.fasta(s));
      try {
        final AssembleParams params = AssembleParams.builder().kmerSize(5).wordSize(4).stepSize(4).reads(Collections.singletonList(input)).directory(tmpDir).create();
        AssembleTask.assemble(params, NullStreamUtils.getNullPrintStream(), new GraphMapStatistics(null));
      } finally {
        FileHelper.deleteAll(input);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
