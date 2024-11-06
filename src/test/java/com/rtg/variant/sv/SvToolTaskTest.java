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
package com.rtg.variant.sv;

import static com.rtg.sam.SharedSamConstants.OUT_SAM;
import static com.rtg.sam.SharedSamConstants.REF_SEQS;
import static com.rtg.sam.SharedSamConstants.SAM_UNSORTED;
import static com.rtg.util.StringUtils.FS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.SimpleArchive;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SvToolTaskTest extends AbstractNanoTest {

  public void testUnsortedReadStats() throws IOException {
    check(REF_SEQS, SAM_UNSORTED, "RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new String[] {"is not sorted in coordinate order."}, 1);
  }

  public void testMakeParams() throws IOException, InvalidParamsException {
    check(REF_SEQS, SAM_UNSORTED, "RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new String[] {"is not sorted in coordinate order."}, 1);
  }

  protected void check(String refSeq, String sam, String readStats,
                       String[] errorMsgs, int errCode) throws IOException {
    try (TestDirectory tdir = new TestDirectory("svtooltask")) {

      final File output = new File(tdir, "sv_out");
      final File input = FileUtils.createTempDir("testcheck", "sv_in", tdir);
      FileUtils.stringToFile(sam, new File(input, OUT_SAM));
      final File readStatsFile = new File(input, "readstats");
      FileUtils.stringToFile(readStats, readStatsFile);
      final String outn = output.getPath();
      final String inn = input.getPath();
      final File templ = ReaderTestUtils.getDNASubDir(refSeq, tdir);
      final MainResult res = MainResult.run(new SvToolCli(), "-t", templ.getPath(), "-o", outn, "-r", readStatsFile.toString(), inn + FS + OUT_SAM, "-s", "1", "-f", "1", "--simple-signals");
      //System.err.println("err\n" + errStr);
      assertEquals("Error" + res.err(), errCode, res.rc());
      //System.err.println(errStr);
      TestUtils.containsAll(res.err(), errorMsgs);
    }
  }

  public void testNormalize() {
    final double[] n1 = SvToolTask.normalize(new double[] {1.0, 2.0, 3.0});
    assertEquals("[0.298, -0.489, -1.005]", Utils.realFormat(n1, 3));
  }

  public void testNormalize1() {
    final double inc = 1000000.0;
    final double[] n1 = SvToolTask.normalize(new double[] {1.0 + inc, 2.0 + inc, 3.0 + inc});
    assertEquals("[0.298, -0.489, -1.005]", Utils.realFormat(n1, 3));
  }

  public void testNormalize2() {
    final double[] n1 = SvToolTask.normalize(new double[] {1.0, 12.0, 22.0});
    assertEquals("[4.777, -4.777, -9.120]", Utils.realFormat(n1, 3));
  }

  public void testNormalize3() {
    final double[] n1 = SvToolTask.normalize(new double[] {1.0, 50.0, 100.0});
    assertEquals("[21.280, -21.280, -42.995]", Utils.realFormat(n1, 3));
  }

  public void testInteresting() throws Exception {
    try (final TestDirectory tempDir = new TestDirectory("svtool")) {
      final File sam = FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvMappings.sam.gz", new File(tempDir, "mappings.sam.gz"));
      FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvMappings.sam.gz.tbi", new File(tempDir, "mappings.sam.gz.tbi"));
      final File rgstats = FileHelper.resourceToFile("com/rtg/variant/sv/resources/rgstats.txt", new File(tempDir, "rgstats.tsv"));
      final File output = new File(tempDir, "output");
      final File template = new File(tempDir, "template");
      assertTrue(template.mkdir());
      SimpleArchive.unpackArchive(FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvTemplate.dwa", new File(tempDir, "template.dwa")), template);

      final MainResult res = MainResult.run(new SvToolCli(),
        "--Xheterozygous",
        "--simple-signals",
        "-r", rgstats.getPath(),
        "-t", template.getPath(),
        "-o", output.getPath(),
        "-s", "10",
        "--region", "simulatedSequence1:47400-58600",
        "--fine-step", "1",
        sam.getPath());
      assertEquals(res.err(), 0, res.rc());

      final String bayes = StringUtils.grep(FileHelper.gzFileToString(new File(output, "sv_bayesian.tsv.gz")), "^[^#]");
      mNano.check("small_sv_bayes.bed", bayes);
      final String interesting = FileHelper.gzFileToString(new File(output, "sv_interesting.bed.gz"));
      mNano.check("small_sv_interesting.bed", interesting);
      final String simple = StringUtils.grep(FileHelper.gzFileToString(new File(output, "sv_simple.tsv.gz")), "^[^#]");
      mNano.check("small_sv_simple.txt", simple);
    }

  }
}
