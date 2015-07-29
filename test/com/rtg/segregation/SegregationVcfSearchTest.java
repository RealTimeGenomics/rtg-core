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

package com.rtg.segregation;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.GlobalFlags;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.ReferenceGenome;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SegregationVcfSearchTest extends TestCase {

  public void check(final String exp, final String expBed, final String res, int newPen, int xoPen) throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      GlobalFlags.resetAccessedStatus();
      Diagnostic.setLogStream();
      final File template = ReaderTestUtils.getDNADir(">1\nACGT\n>2\nACGT\n>3\nACGT\n>X\nACGTGAACGTGAACGTGAACGTGAAAAAACGTGAA\n>Y\nACAAACGTGAACGTGAACGTGAAAAGAA", new File(dir, "template"));
      FileUtils.stringToFile("version 1\n"
          + "either def diploid linear\n"
          + "male seq X haploid linear Y\n"
          + "male seq Y haploid linear Y\n"
          + "female seq X diploid linear\n"
          + "female seq Y none linear\n"
          + "male dup X:1-16 Y:1-16\n"
          , new File(template, ReferenceGenome.REFERENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File vcf = FileHelper.resourceToFile(res, new File(dir, "vcf.vcf"));
      final File output = new File(dir, "out.txt");
      final File regions = new File(dir, "regions.bed");
      final String[] args = {"--vcf", vcf.getPath(), "--template", template.getPath(), "--father", "Father", "--mother", "Mother", "--output", output.getPath(), "--regions-output", regions.getPath(), "-Z", "--Xnew-penalty", "" + newPen, "--Xxo-penalty", "" + xoPen};
      final int ret = new SegregationVcfSearch().mainInit(args, out.printStream(), err.printStream());
      assertEquals(err.toString(), 0, ret);
      //System.err.println(out.toString());
      //System.err.println(err.toString());
      final String actual = IOUtils.readAll(output);
      assertEquals(exp, actual);
      assertEquals("", err.toString());
      final String actualRegions = IOUtils.readAll(regions);
      assertEquals(expBed, actualRegions.replace('\t', ' '));
    }
  }

  public void test() throws IOException {
    Diagnostic.setLogStream();
    final String exp = ""
        + "NS 1" + LS
        + "ME 1 568954 0_0 0_0 0_2 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "HO 1 568968 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "HO 1 645605 0_0 1_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1" + LS
        + "OK 1 705882 0_1 0_0 0_1 0_0 0_0 0_1 0_1 0_0 0_1 0_0 0_1 0_1 0_1" + LS
        + "HO 1 723798 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1" + LS
        + "BN 1 705882 - 705882 1 fa 10011010111 mo ??????????? 10011010111 ???????????" + LS
        + "NS 2" + LS
        + "HO 2 723891 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1" + LS
        + "OK 2 725104 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1" + LS
        + "AH 2 725105 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1 0_1" + LS
        + "OK 2 725106 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1" + LS
        + "BN 2 725104 - 725106 2 fa ??0????0??? mo ??0????0??? ??0????0??? ??0????0???" + LS
        + "NS 3" + LS
        + "OK 3 736289 0_0 0_1 0_0 0_0 0_1 0_1 0_0 0_0 0_1 0_1 0_1 0_0 0_1" + LS
        + "OK 3 740894 0_1 0_0 0_0 0_1 0_1 0_0 0_0 0_1 0_0 0_1 0_0 0_0 0_0" + LS
        + "OK 3 752566 1_1 0_1 1_1 1_1 0_1 0_1 1_1 1_1 0_1 0_1 0_1 1_1 0_1" + LS
        + "OK 3 752567 1_1 0_1 1_1 1_1 1_1 1_1 1_1 1_1 0_1 0_1 0_1 1_1 0_1" + LS
        + "OK 3 752568 1_1 0_1 1_1 1_1 0_1 0_1 1_1 1_1 0_1 0_1 0_1 1_1 0_1" + LS
        + "BN 3 736289 - 752566 3 fa 01100101000 mo 00110011101 01100101000 00110011101" + LS
        + "BE 3 752567 - 752567 1 fa ??????????? mo 11111100010 01100101000 00110011101" + LS
        + "BL 3 752568 - 752568 1 fa ??????????? mo 11001100010 01100101000 00110011101" + LS
        ;
    final String expBed = "3 0 198022430 01100101000 00110011101 N 4 5" + LS;
    check(exp, expBed, "com/rtg/segregation/resources/test.vcf", Search.DEFAULT_NEW_PENALTY, Search.DEFAULT_XO_PENALTY);
  }

  public void testCrossover() throws IOException {
    final String exp = "NS 1" + LS
        + "OK 1 568954 0_2 0_0 0_2 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK 1 568968 0_0 0_1 0_0 0_1 0_1 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK 1 645605 0_1 1_1 0_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1 1_1" + LS
        + "OK 1 705882 0_0 0_1 0_0 0_1 0_1 0_1 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK 1 723798 0_0 0_1 0_0 0_1 0_1 0_1 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK 1 723891 0_2 0_0 0_2 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "BN 1 568954 - 645605 3 fa 10000000000 mo 01100000000 10000000000 01100000000" + LS
        + "BX 1 705882 - 723891 3 fa 10000000000 mo 01110000000 10000000000 01110000000" + LS
        + "XO 1 705882 mo 3 9 8" + LS
        + "NS 2" + LS
        + "OK 2 425105 0_2 0_0 0_2 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK 2 425106 0_0 0_1 0_1 0_1 0_0 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1" + LS
        + "OK 2 425107 0_0 0_1 0_0 0_1 0_0 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1" + LS
        + "OK 2 436289 0_2 0_0 0_0 0_2 0_2 0_2 0_2 0_2 0_2 0_2 0_2 0_2 0_2" + LS
        + "OK 2 740894 0_2 0_0 0_2 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_2" + LS
        + "OK 2 752566 0_0 0_1 0_0 0_1 0_0 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1" + LS
        + "BN 2 425105 - 425106 2 fa 10000000000 mo 11011110111 10000000000 11011110111" + LS
        + "BE 2 425107 - 436289 2 fa 01111111111 mo 01011110111 10000000000 11011110111" + LS
        + "BE 2 740894 - 752566 2 fa 10000000001 mo 01011110111 10000000000 11011110111" + LS
        ;
    final String expBed = ""
        + "1 0 705881 10000000000 01100000000 N 1 2" + LS
        + "1 645605 249250621 10000000000 01110000000 X 1 3" + LS
        + "2 0 243199373 10000000000 11011110111 N 1 2" + LS
        ;
    check(exp, expBed, "com/rtg/segregation/resources/crossover.vcf", 10 , 3);
  }
  public void testXChromosomes() throws IOException {
    final String exp = "NS X" + LS
        + "OK X 2 0_1 0_1 0_1 1_1 0_1 0_0 0_1 0_1 0_1 0_0 1_1 0_1 0_0" + LS
        + "OK X 3 1_1 0_1 0_1 1_1 0_1 0_1 1_1 1_1 0_1 0_1 1_1 1_1 0_1" + LS
        + "OK X 4 0_1 1_1 1_1 1_1 1_1 0_1 0_1 0_1 1_1 0_1 1_1 0_1 0_1" + LS
        + "HO X 5 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK X 6 1_1 0_1 1_1 0_1 1_1 1_1 0_1 0_1 1_1 1_1 0_1 0_1 1_1" + LS
        + "OK X 7 0_1 1_1 1_1 1_1 1_1 0_1 0_1 0_1 1_1 0_1 1_1 0_1 0_1" + LS
        + "OK X 8 0_1 0_0 0_1 0_1 0_1 0_0 0_0 0_0 0_1 0_0 0_1 0_0 0_0" + LS
        + "OK X 9 0_1 0_1 1_1 0_1 1_1 0_1 0_0 0_0 1_1 0_1 0_1 0_0 0_1" + LS
        + "OK X 10 0_1 0_1 1_1 1_1 1_1 0_1 0_0 0_0 1_1 0_1 0_1 0_0 0_1" + LS
        + "OK X 11 1_1 0_1 0_1 0_1 0_1 0_1 1_1 1_1 0_1 0_1 1_1 1_1 0_1" + LS
        + "OK X 12 0_1 0_0 0_1 0_1 0_1 0_0 0_0 0_0 0_1 0_0 0_1 0_0 0_0" + LS
        + "OK X 13 0_1 0_1 1_1 1_1 1_1 0_1 0_0 0_0 1_1 0_1 0_1 0_0 0_1" + LS
        + "HO X 14 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0 0_0" + LS
        + "OK X 15 0_1 0_1 1_1 0_1 0_1 0_1 0_1 0_0 0_1 0_1 0_1 0_1 0_1" + LS
        + "OK X 20 0 0_1 0_1 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "OK X 22 0 0_1 0_1 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "OK X 23 0 0_1 0_1 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "OK X 24 0 0_1 0_1 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "OK X 25 0 0_1 0_1 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "OK X 568924 0 0_1 0_1 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "OK X 568950 0 0_1 0_0 0_0 0_0 0 0 0 0_0 0 0_0 0 0" + LS
        + "BN X 2 - 9 7 fa 11100010100 mo 01001100110 11100010100 01001100110" + LS
        + "BX X 10 - 13 4 fa 11100010100 mo 11110011001 11100010100 00001100110" + LS
        + "XO X 10 mo 1 6 7" + LS
        + "BE X 15 - 568924 7 fa 11111011111 mo 10000000000 11100010100 00001100110" + LS
        + "BE X 568950 - 568950 1 fa 000???0?0?? mo 00000000000 11100010100 00001100110" + LS
        ;
    final String expBed = "X 0 9 11100010100 01001100110 N 5 5" + LS
        + "X 9 243199373 11100010100 00001100110 X 5 4" + LS
        ;
    check(exp, expBed, "com/rtg/segregation/resources/Xchromosome.vcf", 10 , 3);
  }
  public void testPloidyMismatchWarning() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      Diagnostic.setLogStream();
      final File template = ReaderTestUtils.getDNADir(">1\nACGT\n>2\nACGT\n>3\nACGT\n>X\nACGTGAACGTGAACGTGAACGTGAAAAAACGTGAA\n>Y\nACAAACGTGAACGTGAACGTGAAAAGAA", new File(dir, "template"));
      FileUtils.stringToFile("version 1\n"
          + "either def diploid linear\n"
          + "male seq X haploid linear Y\n"
          + "male seq Y haploid linear Y\n"
          + "female seq X diploid linear\n"
          + "female seq Y none linear\n"
          , new File(template, ReferenceGenome.REFERENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File vcf = FileHelper.resourceToFile("com/rtg/segregation/resources/Xchromosome.vcf", new File(dir, "vcf.vcf"));
      final File output = new File(dir, "out.txt");
      final File regions = new File(dir, "regions.bed");
      final String[] args = {"--vcf", vcf.getPath(), "--template", template.getPath(), "--father", "Father", "--mother", "Mother", "--output", output.getPath(), "--regions-output", regions.getPath(), "-Z", "--Xnew-penalty", "10", "--Xxo-penalty", "3"};
      final int ret = new SegregationVcfSearch().mainInit(args, out.printStream(), err.printStream());
      assertEquals(err.toString(), 0, ret);
      TestUtils.containsAll(err.toString(), "There were 14 variants with genotypes that did not match the expected ploidy in chromosome X");
    }
  }
}
