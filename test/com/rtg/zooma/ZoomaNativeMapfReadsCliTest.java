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
package com.rtg.zooma;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.test.RandomDna;
import com.rtg.mode.DnaUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

/**
 */
public class ZoomaNativeMapfReadsCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ZoomaNativeMapfReadsCli();
  }

  public void test() throws IOException {
    if (NativeZooma.isEnabled()) {
      try (TestDirectory dir = new TestDirectory()) {
        final String reference = RandomDna.random(1000);
        final File referenceFile = new File(dir, "reference.fasta");
        FileUtils.stringToFile(">reference\n" + reference + "\n", referenceFile);
        final File indexFile = new File(dir, "reference.bin");
        final MemoryPrintStream mpout = new MemoryPrintStream();
        final MemoryPrintStream mperr = new MemoryPrintStream();
        new ZoomaNativeBuildIndexCli().mainInit(new String[] {"-i", referenceFile.toString(), "-o", indexFile.toString()}, mpout.outputStream(), mperr.printStream());

        // Simple PE mapf
        final File leftFile = new File(dir, "reads_l.fastq");
        FileUtils.stringToFile(""
            + ZoomaNativeMapReadsCliTest.fastqize("read0", reference.substring(100, 200))
            + ZoomaNativeMapReadsCliTest.fastqize("read1", DnaUtils.reverseComplement(reference.substring(350, 450))), leftFile);
        final File rightFile = new File(dir, "reads_r.fastq");
        FileUtils.stringToFile(""
            + ZoomaNativeMapReadsCliTest.fastqize("read0", DnaUtils.reverseComplement(reference.substring(300, 400)))
            + ZoomaNativeMapReadsCliTest.fastqize("read1", reference.substring(150, 250)), rightFile);
        final File outFile = new File(dir, "mapf_pe_out");
        checkMainInitOk("-i", indexFile.toString(), "-o", outFile.toString(), "-l", leftFile.toString(), "-r", rightFile.toString(), "-Q", "-m", "200", "-M", "400");

        // Simple SE mapf
        final File outFile2 = new File(dir, "mapf_se_out");
        checkMainInitOk("-i", indexFile.toString(), "-o", outFile2.toString(), "-l", leftFile.toString(), "-Q", "-m", "200", "-M", "400");
      }
    } else {
      //fail(); //
      System.err.println("skipping zmapf test (no native lib)");
    }
  }
}
