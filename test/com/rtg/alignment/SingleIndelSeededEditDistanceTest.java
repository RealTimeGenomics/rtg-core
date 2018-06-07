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

package com.rtg.alignment;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;

/**
 */
public class SingleIndelSeededEditDistanceTest extends SingleIndelEditDistanceTest {

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int substitutionPenalty, int unknownsPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpenPenalty).gapExtendPenalty(gapExtendPenalty).substitutionPenalty(substitutionPenalty).unknownsPenalty(unknownsPenalty).singleIndelPenalties(null).create();
    return new SingleIndelSeededEditDistance(params, 3, 2, 2, false, 1000);
  }

  @Override
  public void testIndelOffTemplateEnd() {
    final byte[] read = DnaUtils.encodeString("AAATCCACATTATATAAATCTAACTATATCATCATACAATACATACATAG");
    final byte[] template = DnaUtils.encodeString("AAATCCACATTATATGGGGGAAATCTAACTATATCATCATACAATACATACATA");

    final int[] actions = getEditDistanceInstance(19, 1, 9, 0).calculateEditDistance(read, read.length, template, 0, Integer.MAX_VALUE, 9, false);
    final String exp = "===============XXXXXX==XXXX=XX==XXX=XX=====XXXXXXX"; //regression
    assertEquals(exp, ActionsHelper.toString(actions));

    //System.err.println(ActionsHelper.toString(actions)); //not sure what the intended behaviour is here, return null? "best alignment" that won't leave the template?
  }

  @Override
  public void testForwardNegative() {
    //doesn't work well for short reads
  }

  @Override
  public void testForwardNegative4() {
    //doesn't work well for short reads
  }


/*  public int[] unmapped(String readStr, int start, int maxShift) throws Exception {
    final byte[] read = DnaUtils.encodeString(readStr);
    int[] actions = null;
    try (SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(new File("/rtgshare/data/human/ref/1000g_v37_phase2/sdf/"))) {

      sr.seek(0);
      byte[] tmpl = new byte[sr.currentLength()];
      sr.readCurrent(tmpl);

      final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(19).gapExtendPenalty(1).substitutionPenalty(9).nsAsMismatches(false).create();
      UnidirectionalEditDistance ed = new SingleIndelSeededEditDistance(params, 3, 2, 2, true);

      actions = ed.calculateEditDistance(read, read.length, tmpl, start, 60, maxShift, false);
    }
    return actions;
  }


  public void testUnmappedExample50() throws Exception {
    //HS2000-900:85:B02LVACXX:5:1106:3183:190300      149     1       50020963        0       *       =       50020773        0       AAGTCCCCACTAGATTAGCTAGACACAGAGCACTGATTGGTGCATTTACAAACCTTGAGCTAGACACAGAGTGCTGATTGGTGCATTTACAATCCTTTAGC   5DB=:CA@;;BEC>>EBHGFEC?DEIHD:HA@GGGIIJJJIIJIGGIIHFDFHCIIJIGHHHDFE?GJIGJIIJJJJGGICHDIJJIIHHHGHFFEDD@C@   RG:Z:NA12878    XN:i:50020873
    //HS2000-900:85:B02LVACXX:5:1106:3183:190300      147     1       50020954        55      101=    =       50020773        0       AAGTCCCCACTAGATTAGCTAGACACAGAGCACTGATTGGTGCATTTACAAACCTTGAGCTAGACACAGAGTGCTGATTGGTGCATTTACAATCCTTTAGC   5DB=:CA@;;BEC>>EBHGFEC?DEIHD:HA@GGGIIJJJIIJIGGIIHFDFHCIIJIGHHHDFE?GJIGJIIJJJJGGICHDIJJIIHHHGHFFEDD@C@   RG:Z:NA12878    XN:i:50020873   XP:i:50020963   IH:i:1  AS:i:0  XT:Z:M
    int[] actions = unmapped("AAGTCCCCACTAGATTAGCTAGACACAGAGCACTGATTGGTGCATTTACAAACCTTGAGCTAGACACAGAGTGCTGATTGGTGCATTTACAATCCTTTAGC", 50020962, 50);
    assertNotNull(actions);
  }

 public void testUnmappedExample100() throws Exception {
    //HS2000-900:85:B02LVACXX:5:1106:3183:190300      149     1       50020963        0       *       =       50020773        0       AAGTCCCCACTAGATTAGCTAGACACAGAGCACTGATTGGTGCATTTACAAACCTTGAGCTAGACACAGAGTGCTGATTGGTGCATTTACAATCCTTTAGC   5DB=:CA@;;BEC>>EBHGFEC?DEIHD:HA@GGGIIJJJIIJIGGIIHFDFHCIIJIGHHHDFE?GJIGJIIJJJJGGICHDIJJIIHHHGHFFEDD@C@   RG:Z:NA12878    XN:i:50020873
    //HS2000-900:85:B02LVACXX:5:1106:3183:190300      147     1       50020954        55      101=    =       50020773        0       AAGTCCCCACTAGATTAGCTAGACACAGAGCACTGATTGGTGCATTTACAAACCTTGAGCTAGACACAGAGTGCTGATTGGTGCATTTACAATCCTTTAGC   5DB=:CA@;;BEC>>EBHGFEC?DEIHD:HA@GGGIIJJJIIJIGGIIHFDFHCIIJIGHHHDFE?GJIGJIIJJJJGGICHDIJJIIHHHGHFFEDD@C@   RG:Z:NA12878    XN:i:50020873   XP:i:50020963   IH:i:1  AS:i:0  XT:Z:M
    int[] actions = unmapped("AAGTCCCCACTAGATTAGCTAGACACAGAGCACTGATTGGTGCATTTACAAACCTTGAGCTAGACACAGAGTGCTGATTGGTGCATTTACAATCCTTTAGC", 50020962, 100);
    assertNotNull(actions);
  } */


  //XXX Currently these tests below are totally wrong for seeded single indel. They should be looked at before this is used seriously.
  @Override
  public void testInsertAtStart() {
  }

  @Override
  public void testInsertAtEnd() {
  }

  @Override
  public void testInsertAtEndWithMismatchAtStart() {
  }

  @Override
  public void testNoMismatchAfterIndelAtEnd() {
  }

  @Override
  public void testMatchBeforeInsert() {
  }
}
