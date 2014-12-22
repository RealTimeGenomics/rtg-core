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

package com.rtg.simulation;

import com.rtg.sam.SamUtils;
import com.rtg.simulation.SoftClipCigarParser.CigarActionState;

import junit.framework.TestCase;

/**
 */
public class SoftClipCigarParserTest extends TestCase {

  public void testNoSoftClip() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();
    dcp.setState(1, 2, "7=", "7="); //, DnaUtils.encodeString("TTCCCAT"));
    assertEquals(0, dcp.parse());
  }

  public void testSingleSoftClip() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     *  TTCCCAT
     *   TCCCAT
     * TTTCCCATGGTTTTTC
     */
    dcp.setState(1, 2, "7=", "1S6="); //, DnaUtils.encodeString("TTCCCAT"));
    assertEquals(0, dcp.parse());
    /*
     * end gen tpos=2, rpos=1
     * end ali tpos=2, rpos=1
     */
  }

  public void testSingleDel() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * A TCCCATG
     *    CCCATG
     * ATTCCCATGGTTTTTC
     */
    dcp.setState(0, 3, "1=1D7=", "2S6="); //, DnaUtils.encodeString("ATCCCATG"));
    assertEquals(-1, dcp.parse());
    /*
     * end gen tpos=3, rpos=2
     * end ali tpos=3, rpos=2
     */
  }
  public void testSingleDel2() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * A  TCCCATG
     *     CCCATG
     * ATTTCCCATGGTTTTTC
     */
    dcp.setState(0, 4, "1=2D7=", "2S6="); //, DnaUtils.encodeString("ATCCCATG"));
    assertEquals(-2, dcp.parse());
    /*
     * end gen tpos=3, rpos=2
     * end ali tpos=3, rpos=2
     */
  }

  public void testSingleDelEnd() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TT CCCATG
     *    CCCATG
     * TTTCCCATGGTTTTTC
     */
    dcp.setState(0, 3, "2=1D5=", "2S5="); //, DnaUtils.encodeString("TTCCCATG"));
    assertEquals(-1, dcp.parse());
    /*
     * end gen tpos=3, rpos=2
     * end ali tpos=3, rpos=2
     */
  }
  public void testSingleDelEnd2() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TT  CCCATG
     *     CCCATG
     * TTTTCCCATGGTTTTTC
     */
    dcp.setState(0, 4, "2=2D5=", "2S5="); //, DnaUtils.encodeString("TTCCCATG"));
    assertEquals(-2, dcp.parse());
    /*
     * end gen tpos=4, rpos=2
     * end ali tpos=4, rpos=2
     */
  }

  public void testRetardedSoftClipInMiddleOfDel() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TT  CCATG
     *     CCATG
     * TTTCCCATGGTTTTTC
     */
    dcp.setState(0, 3, "2=2D5=", "2S1D5="); //, DnaUtils.encodeString("TTCCCATG"));
    assertEquals(-1, dcp.parse());
    /*
     * end gen tpos=3, rpos=2
     * end ali tpos=3, rpos=2
     */
  }

  public void testDelAndAmbigDel() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     *            1
     * 0     5    0
     * TTCC ATGGTTTT C
     *           TTTTC
     * TTCCCATGGTTTTTC
     * 0    5    1
     *           0
     */
    dcp.setState(0, 10, "4=1D8=1D1=", "8S5="); //, DnaUtils.encodeString("TTCCATGGTTTTC"));
    assertEquals(-2, dcp.parse());

    /*
     * end gen tpos=9, rpos=8
     * end ali tpos=10, rpos=8
     */
  }

  public void testDelAndComplexAmbigDel() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     *            1
     * 0     5    0
     * TTCC ATGGATAT  C
     *            ATATC
     * TTCCCATGGATATATC
     * 0    5    1
     *           0
     */
    dcp.setState(0, 11, "4=1D8=2D1=", "8S5="); //, DnaUtils.encodeString("TTCCATGGATATC"));
    assertEquals(-3, dcp.parse());

    /*
     * end gen tpos=9, rpos=8
     * end ali tpos=10, rpos=8
     */
  }

  public void testIns() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * 0   4
     * ATGTCCCATG
     *     CCCATG
     * AT TCCCATGGTTTTTC
     * 0   3
     */
    dcp.setState(0, 3, "2=1I7=", "4S6="); //, DnaUtils.encodeString("ATGTCCCATG"));
    assertEquals(1, dcp.parse());

    /*
     * end gen tpos=3, rpos=4
     * end ali tpos=3, rpos=4
     */
  }

  public void testIns1() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * ATTTCCCATG
     *     CCCATG
     * ATT CCCATGGTTTTTC
     */
    dcp.setState(0, 3, "3=1I6=", "4S6="); //, DnaUtils.encodeString("ATTTCCCATG"));
    assertEquals(1, dcp.parse());

    /*
     * end gen tpos=3, rpos=4
     * end ali tpos=3, rpos=4
     */
  }

  public void testInsertAmbig() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TTCCCCATG
     *   CCC ATG
     * TTCCC ATGGTTTTTC
     */
    dcp.setState(0, 2, "5=1I3=", "2S6=");
    assertEquals(0, dcp.parse());
  }
  public void testInsertAmbig1() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TTCCCCATG
     *   CC CATG
     * TTCC CATGGTTTTTC
     */
    dcp.setState(0, 2, "4=1I4=", "2S6=");
    assertEquals(0, dcp.parse());
  }
  public void testInsertAmbig2() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TTCCCCATG
     *   C CCATG
     * TTC CCATGGTTTTTC
     */
    dcp.setState(0, 2, "3=1I5=", "2S6=");
    assertEquals(0, dcp.parse());
  }
  public void testInsertAmbig3() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TTCCCCATG
     *    CCCATG
     * TT CCCATGGTTTTTC
     */
    dcp.setState(0, 2, "2=1I6=", "2S6=");
    assertEquals(0, dcp.parse());
  }

  public void testInsertOnBoth() {
    final SoftClipCigarParser dcp = new SoftClipCigarParser();

    /*
     * TTCCCCAATG
     *   CCC AATG
     * TTCCC  ATGGTTTTTC
     */
    dcp.setState(0, 2, "5=2I3=", "2S3=1I3=");
    assertEquals(0, dcp.parse());
  }

  public void testActionsHandler() {
    final CigarActionState cas = new CigarActionState();

    cas.init("1S1X2=1I2=1D1N1=", 3);

    assertEquals(0, cas.getCurrentActionLength());
    assertEquals(0, cas.getCurrentActionRemaining());
    assertEquals(0, cas.getReadPos());
    assertEquals(3, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_SOFT_CLIP, cas.getAction());
    assertEquals(1, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(0, cas.getReadPos());
    assertEquals(3, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_MISMATCH, cas.getAction());
    assertEquals(1, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(1, cas.getReadPos());
    assertEquals(3, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_SAME, cas.getAction());
    assertEquals(2, cas.getCurrentActionLength());
    assertEquals(2, cas.getCurrentActionRemaining());
    assertEquals(2, cas.getReadPos());
    assertEquals(4, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_SAME, cas.getAction());
    assertEquals(2, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(3, cas.getReadPos());
    assertEquals(5, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_INSERTION_INTO_REF, cas.getAction());
    assertEquals(1, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(4, cas.getReadPos());
    assertEquals(6, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_SAME, cas.getAction());
    assertEquals(2, cas.getCurrentActionLength());
    assertEquals(2, cas.getCurrentActionRemaining());
    assertEquals(5, cas.getReadPos());
    assertEquals(6, cas.getTemplatePos());
    assertEquals(SamUtils.CIGAR_SAME, cas.getAction());
    assertEquals(2, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(6, cas.getReadPos());
    assertEquals(7, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_DELETION_FROM_REF, cas.getAction());
    assertEquals(1, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(7, cas.getReadPos());
    assertEquals(8, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_GAP_IN_READ, cas.getAction());
    assertEquals(1, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(7, cas.getReadPos());
    assertEquals(9, cas.getTemplatePos());

    assertEquals(SamUtils.CIGAR_SAME, cas.getAction());
    assertEquals(1, cas.getCurrentActionLength());
    assertEquals(1, cas.getCurrentActionRemaining());
    assertEquals(7, cas.getReadPos());
    assertEquals(10, cas.getTemplatePos());

    assertEquals((char) -1, cas.getAction());
  }

  public void testActionsHandlerBroken() {
    final CigarActionState cas = new CigarActionState();

    cas.init("5", 3);

    assertEquals((char) -1, cas.getAction());
  }
}
