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
package com.rtg.index.hash.ngs.instances;


import java.io.IOException;

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.util.StringUtils;


/**
 */
public class SplitL36w18s2e1Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new SplitL36w18s2e1(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  private static final String PREFIX = ""
    + "Split l=36 w=18 no indels" + StringUtils.LS
    + "number windows=6" + StringUtils.LS
    + "set name=1234 length=23" + StringUtils.LS
    ;
  private static final String SUFFIX = "done" + StringUtils.LS;

  private static final String EXPECTED_1 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000011:00110011:00110010:10101010:10101010" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001001:10011001:10011001:01010101:01010101" + StringUtils.LS
    + "readCall index=2 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001100:11001100:11001110:10101010:10101010" + StringUtils.LS
    + "readCall index=3 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001001:10011010:10101001:10011000:10101010" + StringUtils.LS
    + "readCall index=4 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001100:11001101:01010100:11001101:01010101" + StringUtils.LS
    + "readCall index=5 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001100:11001101:01010101:10011000:10101010" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000011:00110011:00110010:10101010:10101010" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001001:10011001:10011001:01010101:01010101" + StringUtils.LS
    + "templateCall position=7 index=2" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001100:11001100:11001110:10101010:10101010" + StringUtils.LS
    + "templateCall position=7 index=3" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001001:10011010:10101001:10011000:10101010" + StringUtils.LS
    + "templateCall position=7 index=4" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001100:11001101:01010100:11001101:01010101" + StringUtils.LS
    + "templateCall position=7 index=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001100:11001101:01010101:10011000:10101010" + StringUtils.LS
    + SUFFIX
    ;
  public void test1() throws IOException {
    final String dna = "TGCATGCATGCATGCATGCATGCATGCATGCATGCA";
    assertEquals(36, dna.length());
    check(dna, dna, EXPECTED_1);
  }

  private static final String EXPECTED_2 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=2 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=3 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=4 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=5 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=2" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=3" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=4" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + StringUtils.LS
    + SUFFIX
    ;
  public void test2() throws IOException {
    final String dna = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    assertEquals(36, dna.length());
    check(dna, dna, EXPECTED_2);
  }

  private static final String EXPECTED_3 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + StringUtils.LS
    + "readCall index=2 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + StringUtils.LS
    + "readCall index=3 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "readCall index=4 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "readCall index=5 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + StringUtils.LS
    + "templateCall position=7 index=2" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + StringUtils.LS
    + "templateCall position=7 index=3" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=4" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + SUFFIX
    ;
  public void test3() throws IOException {
    final String dna = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    assertEquals(36, dna.length());
    check(dna, dna, EXPECTED_3);
  }

  private static final String EXPECTED_4 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000001:11111111" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111110:00000000" + StringUtils.LS
    + "readCall index=2 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000001:11111111" + StringUtils.LS
    + "readCall index=3 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=4 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000011:11111110:00000000" + StringUtils.LS
    + "readCall index=5 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000011:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000001:11111111" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=2" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000001:11111111" + StringUtils.LS
    + "templateCall position=7 index=3" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=4" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000011:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000011:11111111:11111111" + StringUtils.LS
    + SUFFIX
    ;
  public void test4() throws IOException {
    final String dna = "AAAAAAAAACCCCCCCCCGGGGGGGGGTTTTTTTTT";
    assertEquals(36, dna.length());
    check(dna, dna, EXPECTED_4);
  }

  private static final String EXPECTED_5 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000001:11111111" + StringUtils.LS
    + "readCall index=2 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111110:00000000" + StringUtils.LS
    + "readCall index=3 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000000:00000000:00000000" + StringUtils.LS
    + "readCall index=4 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111111:11111111" + StringUtils.LS
    + "readCall index=5 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111100:00000000:00000000" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000001:11111111" + StringUtils.LS
    + "templateCall position=7 index=2" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111110:00000000" + StringUtils.LS
    + "templateCall position=7 index=3" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000000:00000000:00000000" + StringUtils.LS
    + "templateCall position=7 index=4" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111111:11111111:11111111" + StringUtils.LS
    + "templateCall position=7 index=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000111:11111100:00000000:00000000" + StringUtils.LS
    + SUFFIX
    ;
  public void test5() throws IOException {
    final String dna = "CCCCCCCCCGGGGGGGGGTTTTTTTTTAAAAAAAAA";
    assertEquals(36, dna.length());
    check(dna, dna, EXPECTED_5);
  }

  public void testEmpty() throws IOException {
    check("", "", PREFIX);
  }

  /**
   * Check that all 0, 1, 2 substitutions on the string are found.
   */
  public void testAllSubstitutions() throws IOException {
    final String str = "acgacgtgacacccgtacgtaccccgtgacacccgt";
    assertEquals(36, str.length());
    final Substitute sub = new Substitute(str, SplitL36w18s2e1.FACTORY, true);
    sub.substituteProtected(2);
  }

  public void testFactory() {
    assertEquals(36, SplitL36w18s2e1.FACTORY.hashBits());
    assertEquals(36, SplitL36w18s2e1.FACTORY.windowBits());
    assertEquals(6, SplitL36w18s2e1.FACTORY.numberWindows());
    final NgsHashFunction hf = SplitL36w18s2e1.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof SplitL36w18s2e1);
    assertEquals(hf.numberWindows(), SplitL36w18s2e1.FACTORY.numberWindows());
  }
}

