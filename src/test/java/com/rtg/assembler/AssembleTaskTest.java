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
