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
package com.rtg.metagenomics.metasnp;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.variant.util.arithmetic.SimplePossibility;

public class MetaSnpCliTest extends AbstractCliTest {

  private static final Pattern COMPILE = Pattern.compile(" ");

  @Override
  protected AbstractCli getCli() {
    return new MetaSnpCli();
  }
  public void testValidator() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File inputFile = new File(dir, "input");
      final String input = inputFile.getPath();
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-s", "2", input, "-o", new File(dir, "output").getPath()), "You must provide a value for -m");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-m", "15", "-s", "2", input, "-o", new File(dir, "output").getPath()), "The file '" + input + "' doesn't exist");
      assertTrue(inputFile.createNewFile());
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr("-m", "15", "-s", "1", input, "-o", new File(dir, "output").getPath()), "It makes no sense to run this with less than 2 strains");
      checkHandleFlags("-m", "15", "-s", "2", input, "-o", new File(dir, "output").getPath());
    }
  }
  static final String EXPECTED_VCF = ""
                                 + "##INFO=<ID=LIKE,Number=.,Type=Float,Description=\"phred scaled likelihood of genotype assignments\">\n"
                                 + "##INFO=<ID=SYNDROME,Number=.,Type=String,Description=\"packed representation of strain assignment\">\n"
                                 + "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Strain0 Strain1\n"
                                 + "seq 2 . C A . . LIKE=0.414;SYNDROME=10 GT 1 0\n"
                                 + "seq 9 . G T . . LIKE=0.792;SYNDROME=01 GT 0 1\n"
                                 + "seq 10 . T A . . LIKE=1.139;SYNDROME=10 GT 1 0\n"
                                 + "seq 16 . A C . . LIKE=1.461;SYNDROME=11 GT 1 1\n"
      ;
  static final String[] LINES =  {
  "seq 2 C 1,4 4,1 0,0 0,0 0.0,0.0 0.0,0.0 0.0,0.0 0.0,0.0"
      , "seq 9 G 3,4 4,1 0,0 0,0 0.0,0.0 0.0,0.0 0.0,0.0 0.0,0.0"
      , "seq 10 T 3,3 3,3 0,0 0,0 0.0,0.0 0.0,0.0 0.0,0.0 0.0,0.0"
      , "seq 16 A 8,4 4,1 0,0 0,0 0.0,0.0 0.0,0.0 0.0,0.0 0.0,0.0"
  };

  public void testVcfWriting() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File out = new File(dir, "out.vcf");
      final List<Integer> ref = Arrays.asList(1, 2, 3, 0);
      final List<MetaSnpLine> lines = getLines(LINES);
      final List<AlphaScore> assignments = Arrays.asList(new AlphaScore(0.1, 0.1, 0, 1), new AlphaScore(0.2, 0.2, 2, 3), new AlphaScore(0.3, 0.3, 0, 3), new AlphaScore(0.4, 0.4, 1, 1));
      final EmIterate.EmResult res = new EmIterate.EmResult(new double[][] {{0.1, 0.9}, {0.4, 0.6}}, assignments);
      MetaSnpCli.writeVcf(ref, lines, res, out, SimplePossibility.SINGLETON);
      final String s = FileUtils.fileToString(out).replaceAll("\t", " ");
      assertTrue(s, s.contains(EXPECTED_VCF));
    }

  }
  static final String EXPECTED_VISUAL = "seq 2 C A 2 0.2 0.2" + StringUtils.LS
                                        + "seq 9 G T 1 0.0 0.6" + StringUtils.LS
                                        + "seq 10 T A 2 0.4 0.2" + StringUtils.LS
                                        + "seq 16 A C 3 0.2 0.2" + StringUtils.LS;

  public void testVisualisation() throws IOException {
    try (final ByteArrayOutputStream out = new ByteArrayOutputStream()) {
      final List<Integer> ref = Arrays.asList(1, 2, 3, 0);
      final List<MetaSnpLine> lines = getLines(LINES);
      final List<AlphaScore> assignments = Arrays.asList(new AlphaScore(0.1, 0.1, 0, 1), new AlphaScore(0.2, 0.2, 2, 3), new AlphaScore(0.3, 0.3, 0, 3), new AlphaScore(0.4, 0.4, 1, 1));
      final EmIterate.EmResult res = new EmIterate.EmResult(new double[][] {{0.1, 0.9}, {0.4, 0.6}}, assignments);
      final List<double[][]> evidence = new ArrayList<>();
      evidence.add(new double[][] {{1, 1, 1, 2}, {1, 0, 0, 4}});
      evidence.add(new double[][] {{1, 2, 2, 0}, {0, 1, 1, 3}});
      evidence.add(new double[][] {{2, 1, 2, 0}, {1, 2, 1, 1}});
      evidence.add(new double[][] {{1, 1, 0, 3}, {2, 1, 1, 1}});
      MetaSnpCli.outputVisualisation(ref, lines, evidence, res, out);
      assertEquals(EXPECTED_VISUAL, out.toString().replaceAll("\t", " "));
    }

  }
  public List<MetaSnpLine> getLines(String ... lines) throws IOException {
    final List<MetaSnpLine> list = new ArrayList<>();
    for (String line : lines) {
      list.add(MetaSnpLine.create(COMPILE.matcher(line).replaceAll("\t"), 2));
    }
    return list;
  }

  public void testXiOut() {
    final ByteArrayOutputStream bas = new ByteArrayOutputStream();
    try (PrintStream out = new PrintStream(bas)) {
      MetaSnpCli.writeXi(new double[][]{{0.5, 0.5}, {0.2, 0.8}}, SimplePossibility.SINGLETON, out, Arrays.asList("foo", "bar"));
    }
    assertEquals("\tfoo\tbar" + StringUtils.LS + "Strain0\t0.500\t0.200" + StringUtils.LS + "Strain1\t0.500\t0.800" + StringUtils.LS, bas.toString());
  }
  public void testXiOutNoSamples() {
    final ByteArrayOutputStream bas = new ByteArrayOutputStream();
    try (PrintStream out = new PrintStream(bas)) {
      MetaSnpCli.writeXi(new double[][]{{0.5, 0.5}, {0.2, 0.8}}, SimplePossibility.SINGLETON, out, null);
    }
    assertEquals("\tSample0\tSample1" + StringUtils.LS + "Strain0\t0.500\t0.200" + StringUtils.LS + "Strain1\t0.500\t0.800" + StringUtils.LS, bas.toString());
  }
}
