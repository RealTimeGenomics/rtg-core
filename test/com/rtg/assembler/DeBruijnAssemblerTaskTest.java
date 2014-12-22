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
import java.util.Arrays;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 *         Date: 11/05/12
 *         Time: 2:27 PM
 */
public class DeBruijnAssemblerTaskTest extends TestCase {
  public void test() throws IOException {
    MemoryPrintStream mps = new MemoryPrintStream();
    MemoryPrintStream diag = new MemoryPrintStream();
    File tmpDir = FileHelper.createTempDirectory();
    try {
      File sequence = ReaderTestUtils.getDNADir(""
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCAATGAATCGAAAAAA" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCAATGAATCGAAAAAAGTGTG" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCAATGAATCGAAAAAAGTGTG" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCCATGAATCGAAAAAA" + StringUtils.LS
          + ">a" + StringUtils.LS + "AAAAAACAACCAAGCCATGAATCGAAAAAA" + StringUtils.LS
      );
      try {
        DeBruijnParams params = DeBruijnParams.builder().inputFiles(Arrays.asList(sequence)).directory(tmpDir).kmerSize(6).create();
        File fileList = new File(tmpDir, "fileList");
        FileUtils.stringToFile(sequence.toString() + StringUtils.LS, fileList);
        Diagnostic.setLogStream(diag.printStream());
        try {
          DeBruijnAssemblerTask task = new DeBruijnAssemblerTask(params, mps.outputStream());
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
          Diagnostic.setLogStream();
        }
      } finally {
        FileHelper.deleteAll(sequence);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
