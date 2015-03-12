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

package com.rtg.simulation.sv;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.VcfReaderTest;

import junit.framework.TestCase;
/**
 */
public class VcfReplayerCliTest extends TestCase {

  public VcfReplayerCliTest(String name) {
    super(name);
  }

  private static final String REF1 = ""
    + ">g1" + LS
    // 1234567890123456789012345678901234567890
    + "aaaaaaaaaaccccccccccggggggggggtttttttttt" + LS
    + ">g2" + LS
    + "atatat" + LS
    ;

  private static final String VCF1 = VcfReaderTest.HEADER0
  + "g1" + TAB + "10" + TAB + "." + TAB + "a" + TAB + "a[g1:21[" + TAB + "29" + TAB + "PASS" + TAB + "X=yy;DP=7" + LS
  + "g1" + TAB + "21" + TAB + "." + TAB + "g" + TAB + "[g1:10[g" + TAB + "29" + TAB + "PASS" + TAB + "X=yy;DP=7" + LS
  ;

  private static final String[] RESULT1 = {
   //1234567890123456789012345678901234567890
    "aaaaaaaaaa" +      "ggggggggggtttttttttt",
    "atatat"
  };

  public void test1() throws IOException {
    checkReplay(REF1, VCF1, RESULT1);
  }

  /**
   * This is equivalent to the first (3-break) example from the VCF specification.
   * It exemplifies all possible orientations of breakends in adjacencies.
   */
  private static final String VCF2 = VcfReaderTest.HEADER0
  + "c2     11  bnd_W  G   G]c17:32]  6    PASS  SVTYPE=BND" + LS
  + "c2     12  bnd_V  T   ]c13:26]T  6    PASS  SVTYPE=BND" + LS
  + "c13    26  bnd_U  C   C[c2:12[   6    PASS  SVTYPE=BND" + LS
  + "c13    27  bnd_X  A   [c17:33[A  6    PASS  SVTYPE=BND" + LS
  + "c17    32  bnd_Y  A   A]c2:11]   6    PASS  SVTYPE=BND" + LS
  + "c17    33  bnd_Z  C   [c13:27[C  6    PASS  SVTYPE=BND" + LS
  ;

  private static final String REF2 = ""
    + ">c2" + LS
    // 12345678901234567890
    + "aaaaaaaaaaggggggggga" + LS
    + ">c13" + LS
    // 1234567890123456789012345678901234567890
    + "acacacacacacacacacacacacacatatatatatnnnn" + LS
    + ">c17" + LS
    // 1234567890123456789012345678901234567890
    + "aaaaatttttccccctttttgggggtttttaaccaaccaa" + LS
    ;



  private static final String[] RESULT2 = {
    "aaaaaaaaaag" + reverse("aaaaatttttccccctttttgggggtttttaa"),
    "acacacacacacacacacacacacac" + "gggggggga",
    reverse("atatatatatnnnn") + "ccaaccaa",   // doesn't matter if some outputs are reversed
  };



  public void test2() throws IOException {
    checkReplay(REF2, VCF2.replaceAll("  *", TAB), RESULT2);
  }

  private void checkReplay(String refStr, String vcfStr, String[] resultStr) throws IOException {
    checkReplay(refStr, vcfStr, resultStr, null);
  }

  private void checkReplay(String refStr, String vcfStr, String[] resultStr, String[] resultStrB) throws IOException {
    Diagnostic.setLogStream();
    final File temp = FileUtils.createTempDir("vcfreplay", "test");
    try {
      final File ref = new File(temp, "ref");
      ReaderTestUtils.getDNADir(refStr, ref);
      final File vcf = new File(temp, "vcf");
      FileUtils.stringToFile(vcfStr, vcf);
      final File out = new File(temp, "out");
      final CFlags flags = new CFlags();
      VcfReplayerCli.getFlags(flags);
      assertTrue(flags.setFlags("-t", ref.getPath(),
        "-v", vcf.getPath(),
        "-o", out.getPath()));

      new VcfReplayerCli().process(flags);
      final File left = new File(out, "left");
      assertTrue(left.exists());
      //   printAll(r);
      try (SequencesReader r = SequencesReaderFactory.createMemorySequencesReader(left, true, LongRange.NONE)) {
        assertEquals(resultStr.length, r.numberSequences());
        for (int seq = 0; seq < resultStr.length; seq++) {
          final byte[] b = r.read(seq);
          assertEquals(resultStr[seq].toUpperCase(Locale.getDefault()), DnaUtils.bytesToSequenceIncCG(b));
        }
      }
      final File right = new File(out, "right");
      assertTrue(right.exists());
      if (resultStrB != null) {
        //   printAll(r);
        try (SequencesReader r2 = SequencesReaderFactory.createMemorySequencesReader(right, true, LongRange.NONE)) {
          assertEquals(resultStrB.length, r2.numberSequences());
          for (int seq = 0; seq < resultStrB.length; seq++) {
            final byte[] b = r2.read(seq);
            assertEquals(resultStrB[seq].toUpperCase(Locale.getDefault()), DnaUtils.bytesToSequenceIncCG(b));
          }
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

//  private void printAll(SequencesReader r) {
//    while (r.nextSequence()) {
//      final byte[] b = new byte[r.currentLength()];
//      r.readCurrent(b);
//      System.err.println(DnaUtils.bytesToSequenceIncCG(b));
//    }
//    r.seek(0);
//  }

  private static String reverse(String s) {
    return DnaUtils.reverseComplement(s);
  }

  private static final String REF3 = ">c2" + LS
  // 12345678901234567890123456789012345678901234567890
  + "aaaaaaaaaaaaaaaaaaaaaaaaaaatgggggggggggggggggggggg" + LS
  + ">c13" + LS
  + "cccccccccccccccccgggccccgccccttttttttttttttttttttt" + LS;

  private static final String VCF3 = VcfReaderTest.HEADER0
  + "c2      28  bnd_U  G   ]c13:17]AGTNNNNNCAT     6    PASS  SVTYPE=BND;MATEID=bnd_V" + LS
  + "c13     17  bnd_V  C   CAGTNNNNNCA[c2:28[      6    PASS  SVTYPE=BND;MATEID=bnd_U" + LS;

  private static final String[] RESULT3 = {
    "aaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "cccccccccccccccc" + "CAGTNNNNNCA" + "tgggggggggggggggggggggg",
    reverse("gggccccgccccttttttttttttttttttttt"),   // doesn't matter if some outputs are reversed
  };

  public void test3() throws IOException {
    checkReplay(REF3, VCF3.replaceAll("  *", TAB), RESULT3);
  }


  private static final String REF4 = ">c1" + LS
  // 12345678901234567890123456789012345678901234567890
  + "aaaaaaaaaaagaaaaaaaaaaaaaaatgggggggggggggggggggggg" + LS
  + ">c2" + LS
  + "caaaaaaaaaaaaaaaaaagtaaaaaatgggggggggggacggggggggt" + LS
  ;

  // Delete the middle section of c1 and invert the middle section of c2.
  private static final String VCF4 = VcfReaderTest.HEADER0
  + "c1      12  bnd_M  G   G[c1:28[   6    PASS  SVTYPE=BND;MATEID=bnd_L" + LS
  + "c1      28  bnd_L  T   [c1:12[T   6    PASS  SVTYPE=BND;MATEID=bnd_M" + LS
  // now invert c2:21..40
  + "c2      20  bnd_W  G   G]c2:40]    6    PASS  SVTYPE=BND;MATEID=bnd_U;EVENT=INV0" + LS
  + "c2      21  bnd_V  T   [c2:41[T    6    PASS  SVTYPE=BND;MATEID=bnd_X;EVENT=INV0" + LS
  + "c2      40  bnd_U  A   A]c2:20]    6    PASS  SVTYPE=BND;MATEID=bnd_W;EVENT=INV0" + LS
  + "c2      41  bnd_X  C   [c2:21[C    6    PASS  SVTYPE=BND;MATEID=bnd_V;EVENT=INV0" + LS
  ;

  private static final String[] RESULT4 = {
    "aaaaaaaaaaa" + "G" + "tgggggggggggggggggggggg",
    "caaaaaaaaaaaaaaaaaag" + reverse("taaaaaatggggggggggga") + "cggggggggt",
  };

  public void test4() throws IOException {
    checkReplay(REF4, VCF4.replaceAll("  *", TAB), RESULT4);
  }

  private static final String REF5 = ""
    + ">c2" + LS
    // 12345678901234567890
    + "aaaaaaaaaaggggggggga" + LS
    + ">c13" + LS
    // 1234567890123456789012345678901234567890
    + "acacacacacacacacacacacacacatatatatatnnnn" + LS
    + ">c17" + LS
    // 1234567890123456789012345678901234567890
    + "aaaaatttttccccctttttgggggtttttaaccaaccaa" + LS
    ;

  private static final String VCF5 = VcfReaderTest.HEADER0_B
          + "c2      9 bnd_V  A   ]c13:19]A               6    PASS  SVTYPE=BND;MATEID=bnd_U  GT  1/0" + LS //we're being bad and phasing this using unphased terminology
          + "c13     19 bnd_U  A   A[c2:9[,A[c17:20[   6    PASS  SVTYPE=BND;MATEID=bnd_V,bnd_Z  GT 1/2" + LS
          + "c17     20 bnd_Z  T   ]c13:19]T               6    PASS  SVTYPE=BND;MATEID=bnd_U  GT 0/1" + LS
          ;

  private static final String[] RESULT5 = {
    "aaaaaaaa", //remaining c2
    "acacacacacacacacaca" + "aaggggggggga",
    "aaaaatttttccccctttttgggggtttttaaccaaccaa", //remaining c17
    reverse("cacacacatatatatatnnnn"), //remaining c13
  };

  private static final String[] RESULT5_B = {
    "aaaaaaaaaaggggggggga", //remaining c2
    "acacacacacacacacaca" + "tgggggtttttaaccaaccaa",
    "aaaaatttttccccctttt", //remaining c17
    reverse("cacacacatatatatatnnnn"), //remaining c13
  };

  public void test5() throws IOException {
    checkReplay(REF5, VCF5.replaceAll("  *", TAB), RESULT5, RESULT5_B);
  }
}
