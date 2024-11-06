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

  public void test() throws IOException {
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
