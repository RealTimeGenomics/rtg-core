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
import java.util.Arrays;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.mode.DnaUtils;
import com.rtg.sam.SamUtils;
import com.rtg.util.PortableRandom;
import com.rtg.util.SpawnJvm;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.RandomDna;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 */
public class ZoomaNativeMapReadsCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ZoomaNativeMapReadsCli();
  }

  static String fastqize(String name, String seq) {
    return "@" + name + StringUtils.LS + seq + StringUtils.LS + "+" + name + StringUtils.LS + StringUtils.repeat('5', seq.length()) + StringUtils.LS;
  }

  // Inhale a small sam or bam file into a string
  static String samFileToString(File samFile) throws IOException {
    final StringBuilder samText = new StringBuilder();
    try (final SamReader reader = SamUtils.makeSamReader(samFile)) {
      for (final SAMRecord rec : reader) {
        samText.append(rec.getSAMString());
      }
    }
    return samText.toString();
  }

  // called by main() in the spawned JVM
  private void runActualTest() throws IOException {
    setUp();
    try (TestDirectory dir = new TestDirectory()) {
      final PortableRandom rand = new PortableRandom(42);
      final StringBuilder referenceFa = new StringBuilder();
      final StringBuilder leftFq = new StringBuilder();
      final StringBuilder rightFq = new StringBuilder();
      int read = 0;
      final int reflen = 10000;
      final int readlen = 100;
      final int maxinslen = 200;
      final int maxfraglen = readlen * 2 + maxinslen;
      for (int ref = 0; ref < 5; ++ref) {
        final String reference = RandomDna.random(reflen, rand);
        referenceFa.append(">reference").append(ref).append("\n").append(reference).append("\n");

        do {
          final int pos = rand.nextInt(reflen - readlen - maxfraglen);
          final int posr = pos + maxfraglen - readlen - rand.nextInt(readlen);
          String lbases = reference.substring(pos, pos + readlen);
          String rbases = reference.substring(posr, posr + readlen);
          rbases = DnaUtils.reverseComplement(rbases);
          String orientation = "F1R2";
          if (rand.nextBoolean()) {
            orientation = "F2R1";
            final String t = rbases;
            rbases = lbases;
            lbases = t;
          }
          leftFq.append(fastqize("read" + read + "/1/" + orientation, lbases));
          rightFq.append(fastqize("read" + read + "/2/" + orientation, rbases));
        } while (read++ <= ref * 5);
      }

      final File referenceFile = new File(dir, "reference.fasta");
      FileUtils.stringToFile(referenceFa.toString(), referenceFile);

      final File indexFile = new File(dir, "reference.bin");
      final MemoryPrintStream mpout = new MemoryPrintStream();
      final MemoryPrintStream mperr = new MemoryPrintStream();

      new ZoomaNativeBuildIndexCli().mainInit(new String[] {"-i", referenceFile.toString(), "-o", indexFile.toString()}, mpout.outputStream(), mperr.printStream());

      // Simple PE map
      final File leftFile = new File(dir, "reads_l.fastq");
      FileUtils.stringToFile(leftFq.toString(), leftFile);
      final File rightFile = new File(dir, "reads_r.fastq");
      FileUtils.stringToFile(rightFq.toString(), rightFile);
      final File outFile = new File(dir, "map_out");

      final MainResult result = MainResult.run(getCli(), "-t", indexFile.toString(), "-o", outFile.toString(), "-l", leftFile.toString(), "-r", rightFile.toString(), "-Q", "-T", "1", "-m", "200", "-M", "400", "-N");
      assertEquals(result.err(), 0, result.rc());

      assertTrue(new File(outFile, "done").exists());
      // Concatenate output files
      final StringBuilder sb = new StringBuilder();
      final File[] outFiles = FileUtils.listFiles(outFile, f -> f.getName().endsWith(".bam"));
      Arrays.sort(outFiles);
      for (File f : outFiles) {
        //System.err.println(FileUtils.fileToString(f));
        sb.append(TestUtils.stripSAMHeader(samFileToString(f)));
      }
      mNano.check("testzmap.sam", sb.toString(), false);
    }
  }

  // invoked by spawn in test() method
  public static void main(final String[] args) throws IOException {
    new ZoomaNativeMapReadsCliTest().runActualTest();
  }

  public void test() throws IOException, InterruptedException {
    if (NativeZooma.isEnabled()) {
      final Process p = SpawnJvm.spawn(ZoomaNativeMapReadsCliTest.class.getName());
      final String out = FileUtils.streamToString(p.getErrorStream());
      TestUtils.containsAll(out, "sort time", "total runtime", "max nucleotides to hash");
      assertEquals("", FileUtils.streamToString(p.getInputStream()));
      p.getOutputStream().close();
      p.getInputStream().close();
      p.getErrorStream().close();
    } else {
      System.err.println("skipping zmap test (no native lib)");
    }
  }
}
