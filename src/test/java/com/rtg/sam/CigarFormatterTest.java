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
package com.rtg.sam;

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 * Test class for {@link CigarFormatter}
 */
public class CigarFormatterTest extends TestCase {

  public CigarFormatterTest(final String name) {
    super(name);
  }

  // test different CG overlap sizes
  public void testActionsCgNewOverlap() {
    //tattttcgcaggacttattttaatt.....ctcaaacgct
    for (int overlap = 0; overlap <= 4; ++overlap) {
      final String os = "BBBB".substring(4 - overlap);
      final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=====" + os + "====================NNNNN==========", 0, 0), false, 40, false, true);
      assertEquals((25 - overlap) + "=5N10=", cigar);
    }
  }

  public void testActionsCgNewGapLeftArmNoRC() {
    // tattttcgcaggacttattttaatt....ctcaaacgct
    String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNN==========", 0, 0), false, 40, false, true);
    assertEquals("25=4N10=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNNN==========", 0, 0), false, 40, false, true);
    assertEquals("25=5N10=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("===XB=====================NNNNNN==========", 0, 0), false, 39, false, true);
    assertEquals("3=1X20=6N9=1S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=====BBX===================NNNNNNN==========", 0, 0), false, 38, true, true);
    assertEquals("23M7N8M2S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("===D==BBBB====================NNNNNNNN==========", 0, 0), false, 36, false, true);
    assertEquals("3=1D18=8N6=4S", cigar);
    // this is tricky - it should skip the ==I=, but not the S.
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=====BBBB==D=X================NNNNNNNN==========", 0, 0), false, 36, false, true);
    assertEquals("5=1X16=8N6=4S", cigar);
    // another tricky one - it should skip the ===III=
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=====BBBB===III==============NNNNNNNN==========", 0, 0), false, 35, false, true);
    assertEquals("18=8N9=1S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNN==========", -1, 0), false, 400, false, true);
    assertEquals("1S24=4N10=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNN==========", -3, 0), false, 35, false, true);
    assertEquals("3S22=4N9=1S", cigar);
  }

  public void testActionsCgNewGapLeftArmRc() {
    //tattttcgcaggacttattttaatt.....ctcaaacgct
    assertEquals("10=4N25=", CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNN==========", 0, 0), true, 40, false, true));
    assertEquals("10=5N25=", CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNNN==========", 0, 0), true, 40, false, true));
    assertEquals("10=6N20=1X2=1S", CigarFormatter.actionsToCigar(ActionsHelper.build("===XB=====================NNNNNN==========", 0, 0), true, 39, false, true));
    assertEquals("10M7N21M2S", CigarFormatter.actionsToCigar(ActionsHelper.build("=====BBX===================NNNNNNN==========", 0, 0), true, 38, true, true));
    assertEquals("10=8N17=1D4S", CigarFormatter.actionsToCigar(ActionsHelper.build("====D=BBBB====================NNNNNNNN==========", 0, 0), true, 36, false, true));
//    // this is tricky - it should skip the D as well as the three = and one X.
    assertEquals("10=8N16=1X1=4S", CigarFormatter.actionsToCigar(ActionsHelper.build("=====BBBB==D=X================NNNNNNNN==========", 0, 0), true, 36, false, true));
//    // another tricky one - it should skip the two = then two of the Is
    assertEquals("10=8N13=5S", CigarFormatter.actionsToCigar(ActionsHelper.build("=====BBBB===III==============NNNNNNNN==========", 0, 0), true, 31, false, true));
    assertEquals("1S9=4N25=", CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNN==========", -1, 0), true, 400, false, true));
    assertEquals("3S7=4N24=1S", CigarFormatter.actionsToCigar(ActionsHelper.build("=========================NNNN==========", -3, 0), true, 35, false, true));
  }

  public void testActionsCgNewGapRightArmNoRc() {
    //tattttcgca.....ggacttattttaattctcaaacgct
    String cigar;
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNN=========================", 0, 0), false, 40, false, false);
    assertEquals("10=4N25=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNN=========================", 0, 0), false, 40, false, false);
    assertEquals("10=5N25=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNN=====================BX===", 0, 0), false, 39, false, false);
    assertEquals("10=6N20=1X2=1S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNN===================XBB=====", 0, 0), false, 38, true, false);
    assertEquals("10M7N21M2S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN====================BBBB=D====", -4, 0), false, 36, false, false);
    assertEquals("4S6=8N17=1D4=", cigar);
    // this is tricky - it should skip the I as well as the three = and one X.
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN================X=D==BBBB=====", 0, 0), false, 36, false, false);
    assertEquals("10=8N16=1X1=4S", cigar);
    // another tricky one - it should skip the two = then two of the Is
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN===============III==BBBB=====", 0, 0), false, 31, false, false);
    assertEquals("10=8N13=5S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNN=========================", -1, 0), false, 400, false, false);
    assertEquals("1S9=4N25=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNN=========================", -3, 0), false, 35, false, false);
    assertEquals("3S7=4N24=1S", cigar);
  }

  public void testActionsCgNewGapRightArmRc() {
    //tattttcgca.....ggacttattttaattctcaaacgct
    String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNN=========================", 0, 0), true, 40, false, false);
    assertEquals("25=4N10=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNN=========================", 0, 0), true, 40, false, false);
    assertEquals("25=5N10=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNN=====================BX===", 0, 0), true, 39, false, false);
    assertEquals("3=1X20=6N9=1S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNN===================XBB=====", 0, 0), true, 38, true, false);
    assertEquals("23M7N8M2S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN====================BBBB=D====", -4, 0), true, 36, false, false);
    assertEquals("4S1D17=8N10=", cigar);
    // this is tricky - it should skip the D as well as the three = and one S.
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN================X=D==BBBB=====", 0, 0), true, 36, false, false);
    assertEquals("5=1X16=8N6=4S", cigar);
    // another tricky one - it should skip the two = then two of the Is
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN===============III==BBBB=====", 0, 0), true, 31, false, false);
    assertEquals("18=8N5=5S", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNN=========================", -1, 0), true, 400, false, false);
    assertEquals("1S24=4N10=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNN=========================", -3, 0), true, 35, false, false);
    assertEquals("3S22=4N9=1S", cigar);
  }

  public void testCgAdjacentIDLeftArm() {
    //tattttcgcaggacttattttaatt.....ctcaaacgct
    String cigar;
    // deleting these overlaps puts a I next to an D, which is not allowed.
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("====IBBB===D================NNNNNNNN==========", 0, 0), false, 40, false, true);
    // was: assertEquals("4=1I1D16=8N10=", cigar);
    assertEquals("4=1X16=8N10=", cigar);
    // deleting these overlaps puts a I next to an D, which is not allowed.
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("===IIBBB==DDDD==============NNNNNNNN==========", 0, 0), false, 40, false, true);
    // was: assertEquals("3=2I3D14=8N10=", cigar);
    assertEquals("3=1I1X2D14=8N10=", cigar);
    // deleting these overlaps puts a D next to an I, which is not allowed.
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("=====DBBB==XIII==============NNNNNNNN==========", 0, 0), false, 40, false, true);
    // was: assertEquals("5=1D3I14=8N10=", cigar);
    assertEquals("5=1X2I14=8N10=", cigar);
  }

  public void testCgAdjacentIDRightArm() {
    //tattttcgca.....ggacttattttaattctcaaacgct

    //t     tattttcgcaNNNNNNNNggacttattttaattc---tcaaacgct
    //r     ==========NNNNNNNN===============-   ====
    //                                        III==
    //new                                    X II==

    // was: assertEquals("4S6=8N15=1D3I2=", cigar);
    assertEquals("10=8N16=1X2I2=", CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN================D====BBBBIII==", 0, 0), false, 50, false, false));


    //t     tattttcgcaNNNNNNNNggacttattttaa---ttctcaaacgct
    //r     ==========NNNNNNNN=============III====
    //                                       --=====

    // was: assertEquals("4S6=8N15=3I2D5=", cigar);
    assertEquals("10=8N13=2I1X1D5=", CigarFormatter.actionsToCigar(ActionsHelper.build("==========NNNNNNNN=============III====BBBBDD=====", 0, 0), false, 50, false, false));
  }

  public void testActionsSimple() {
    String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("===III==DD==XXXX==", 0, 0), false, 18, false, true);
    assertEquals("3=3I2=2D2=4X2=", cigar);
    cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("===III==DD==XXXX==", 0, 0), false, 18, true, true);
    assertEquals("3M3I2M2D8M", cigar);
  }

  public void testActionsIdentical() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==================", 0, 0), false, 18, false, true);
    assertEquals("18=", cigar);
  }

  public void testActionsInsertFirst() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("D=================", 0, 0), false, 18, false, true);
    assertEquals("1D17=", cigar);
  }

  public void testActionsSoftClipStart() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==================", -3, 0), false, 18, false, true);
    assertEquals("3S15=", cigar);
  }

  public void testActionsSoftClipEnd() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==================", 3, 0), false, 18, false, true);
    assertEquals("15=3S", cigar);
  }

  public void testActionsSoftClipStartSingle() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==================", -1, 0), false, 18, false, true);
    assertEquals("1S17=", cigar);
  }

  public void testActionsSoftClip1EndSingle() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==================", 1, 0), false, 18, false, true);
    assertEquals("17=1S", cigar);
  }

  public void testActionsLargeSoftClippedOverlapStart() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==============================================================================I===", -26, 0), false, 53, false, true);
    assertEquals("26S52=1I1=2S", cigar);
  }
  public void testActionsLargeSoftClippedOverlapEnd() {
    final String cigar = CigarFormatter.actionsToCigar(ActionsHelper.build("==I==D==========================================================================", 0, 0), false, 32, false, true);
    assertEquals("2=1I2=1D27=47S", cigar);
  }

  public void testCigarRefLength() {
    assertEquals(8, SamUtils.cigarRefLength("8M"));
    assertEquals(8, SamUtils.cigarRefLength("4M1I4M"));
    assertEquals(8, SamUtils.cigarRefLength("3M1D4M"));
    assertEquals(8, SamUtils.cigarRefLength("3X1D4="));
    assertEquals(4, SamUtils.cigarRefLength("3M1N4P"));
    assertEquals(8, SamUtils.cigarRefLength("4M1S4M"));
    assertEquals(0, SamUtils.cigarRefLength("13S23I"));
  }

  public void testSoftClipCigarError() {
//    final StringBuilder sb = new StringBuilder();
//    sb.append("9=1X100=1D9=");
//    try {
//      CigarFormatter.softClipCigar(sb, 10, true);
//      fail();
//    } catch (final RuntimeException e) {
//      assertEquals("Cannot soft clip 10 from start of 9=1X100=1D9=", e.getMessage());
//    }
//    try {
//      CigarFormatter.softClipCigar(sb, 10, false);
//      fail();
//    } catch (final RuntimeException e) {
//      assertEquals("Cannot soft clip 10 from end of 9=1X100=1D9=", e.getMessage());
//    }
  }

  public void testSoftClipCigar() {
    int[] actions = ActionsHelper.build("=====================================================================================================", -42, 192);
    assertEquals("42S59=", CigarFormatter.actionsToCigar(actions, false, 200, false, false));
    assertEquals("42S58=1S", CigarFormatter.actionsToCigar(actions, false, 58, false, false));

    actions = ActionsHelper.build("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX", -1, 192);
    assertEquals("1S100X", CigarFormatter.actionsToCigar(actions, false, 200, false, false));
    assertEquals("1S58X42S", CigarFormatter.actionsToCigar(actions, false, 58, false, false));
  }

  public void testInvalidSoftClippingCG() {
    //GAGGGTTAGG.....GTGAGGGTTTGGGTTAGGGTATTAG
    byte[] tmpl = DnaUtils.encodeString("GAGGGTTAGGGTTAGGGTGAGGGTTAGGGTTAGGG".replaceAll(" ", ""));
    int[] actions = ActionsHelper.build("==========NNNNNN=========X==========B=====", 0, 3);
    //AlignmentResult ar = new AlignmentResult(read, actions, tmpl.length, tmpl, false);
    assertEquals("10=6N9=1X9=5S", CigarFormatter.actionsToCigar(actions, false, tmpl.length, false, false));

    //                            GAGGGTTAGG.....GTGAGGGTTTGGGTTAGGGTATTAG
    tmpl = DnaUtils.encodeString("AAAAAAAAAAGAGGGTTAGGGTTAGGGTGAGGGTTAGGGTTAGGG".replaceAll(" ", ""));
    actions = ActionsHelper.build("==========NNNNNN=========X==========B=====", 10, 3);
    //ar = new AlignmentResult(read, actions, tmpl.length, tmpl, false);
    assertEquals("10=6N9=1X9=5S", CigarFormatter.actionsToCigar(actions, false, tmpl.length, false, false));
  }

  public void testMentalRead() {
    final byte[] tmpl = DnaUtils.encodeString("CTCTAACTTCAnTnnGAAGGAGTGGCACTTCCACCTGCCTCAGCTCATGCGTGATATCCAGGTGGGGGCCCAAGATGGTGTCTTGGAGTCTGGGGTAATGCTTGGAGACAGGGA");
    //    byte[] read = DnaUtils.encodeString("----------TGAGCTGAGG.....GAAGTGCCACTCCTTCACTCCTTAT");   //rc of original read, TGAGCTGAGGGAAGTGCCACTCCTTCACTCCTTAT
    //                                                   A TAAG
    //                                                  GAGT  GAAGGAGTGGCACTTC      CCTCAGCTCA
    final int[] actions = ActionsHelper.build("==========NNNNNN=================I=XBBBB=II==", 10, 8);
    assertEquals(44, ActionsHelper.zeroBasedTemplateEndPos(actions));

    assertEquals("2=2I16=6N10=", CigarFormatter.actionsToCigar(actions, true, tmpl.length, false, false));
  }

  public void testRcSoftClipEnd() {
    //CCCAGGAGACGGAGCATCTCCTTCTAATGAGAAACGGTGTAAGGTCTCAGTTAGGATTGAATTTCTTCTTTCTGCACTAAATTGGT
    final byte[] tmpl = DnaUtils.encodeString("GACTAATTTAGTGCAGAAAGAAGAAATTCAATCCTAACTGAGACCTTACACCGTTTCTCATTAGAAGGAGATGCTCCTGTCTCCTGG".replaceAll(" ", ""));
    final int[] actions = ActionsHelper.build("X=========D=========================================================================X==", 1, 4);
    //AlignmentResult ar = new AlignmentResult(read, actions, tmpl.length, tmpl, false);
    assertEquals("2=1X73=1D9=1S", CigarFormatter.actionsToCigar(actions, true, tmpl.length, false, true));
  }

  public void testActionSoftClipStartWithNOOP() {
    final byte[] tmpl = DnaUtils.encodeString("TGCTGGGGAGGACAGGGCAGGGAAGGGCAGGGGGCCAGCAGGAGTGCTGTGGCCGTCCAGACGAGGCCACCCAGACCAGCGGGTCATGAGCTGCTGAACCCTAAAGGCTGGGGTCAGG");
    final int[] actions = ActionsHelper.build("====DDD=================================================================================================", 0, 3);
    ActionsHelper.softClip(actions, true, 4, 3);
    assertEquals("4S97=", CigarFormatter.actionsToCigar(actions, false, tmpl.length, false, true));
  }

  public void testActionSoftClipEndWithNOOP() {
    final byte[] tmpl = DnaUtils.encodeString("TGCTGGGGAGGACAGGGCAGGGAAGGGCAGGGGGCCAGCAGGAGTGCTGTGGCCGTCCAGACGAGGCCACCCAGACCAGCGGGTCATGAGCTGCTGAACCCTAAAGGCTGGGGTCAGG");
    final int[] actions = ActionsHelper.build("=================================================================================================DDD====", 0, 3);
    ActionsHelper.softClip(actions, false, 4, 3);
    assertEquals("97=4S", CigarFormatter.actionsToCigar(actions, false, tmpl.length, false, true));
  }

  public void testZeroLengthCigarNoExplodey() {
    final int[] actions = ActionsHelper.build("", -2, 0);
    final String cigar = CigarFormatter.actionsToCigar(actions, false, 10, false, true);
    assertNotNull(cigar);
    assertEquals(0, cigar.length());
  }

  public void testAlreadySoftClipped() {
//
//    start = 305
//            template length = 372
//    13S55=4S

//    t XXXXXGCCATACCTATCCATCCATCCATGCTATTGGCGTGACCCGACACAGGCGGGACCTGTAAGCCACCGTCAAGTGTATTTCTT
//    r      GGCAGCCCGACACATCCATCCATGCTATTGGCGTGACCCGACACAGGCGGGACCTGTAAGCCACCGTCAAGT
    //       SSSSSSSSSSSSS=======================================================????

    int[] actions = ActionsHelper.build("SSSSSSSSSSSSS===========================================================", 18, 10);
    assertEquals("13S59=", CigarFormatter.actionsToCigar(actions, false, 86, false, true));


    actions = ActionsHelper.build("===========================================================SSSSSSSSSSSSS", 18, 10);
    assertEquals("13S59=", CigarFormatter.actionsToCigar(actions, true, 86, false, true));
  }
}
