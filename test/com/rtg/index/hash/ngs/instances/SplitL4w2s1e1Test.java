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
public class SplitL4w2s1e1Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new SplitL4w2s1e1(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  private static final String PREFIX = ""
    + "Split l=4 w=2 no indels" + StringUtils.LS
    + "number windows=2" + StringUtils.LS
    + "set name=1234 length=23" + StringUtils.LS
    ;
  private static final String SUFFIX = "done" + StringUtils.LS;

  private static final String EXPECTED_1 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000010" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001110" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000010" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001110" + StringUtils.LS
    + SUFFIX
    ;
  public void test1() throws IOException {
    final String dna = "TGCA";
    assertEquals(4, dna.length());
    check(dna, dna, EXPECTED_1);
  }

  private static final String EXPECTED_2 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001111" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001111" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001111" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001111" + StringUtils.LS
    + SUFFIX
    ;
  public void test2() throws IOException {
    final String dna = "TTTT";
    assertEquals(4, dna.length());
    check(dna, dna, EXPECTED_2);
  }

  private static final String EXPECTED_3 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001100" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001100" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001100" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001100" + StringUtils.LS
    + SUFFIX
    ;
  public void test3() throws IOException {
    final String dna = "GGGG";
    assertEquals(4, dna.length());
    check(dna, dna, EXPECTED_3);
  }

  private static final String EXPECTED_4 = ""
    + PREFIX
    + "readCall index=0 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000011" + StringUtils.LS
    + "readCall index=1 id=5" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000011" + StringUtils.LS
    + "templateCall position=7 index=0" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000011" + StringUtils.LS
    + "templateCall position=7 index=1" + StringUtils.LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000011" + StringUtils.LS
    + SUFFIX
    ;
  public void test4() throws IOException {
    final String dna = "CCCC";
    assertEquals(4, dna.length());
    check(dna, dna, EXPECTED_4);
  }

  public void testEmpty() throws IOException {
    check("", "", PREFIX);
  }

  public void testFactory() {
    assertEquals(4, SplitL4w2s1e1.FACTORY.hashBits());
    assertEquals(4, SplitL4w2s1e1.FACTORY.windowBits());
    assertEquals(2, SplitL4w2s1e1.FACTORY.numberWindows());
    final NgsHashFunction hf = SplitL4w2s1e1.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof SplitL4w2s1e1);
    assertEquals(hf.numberWindows(), SplitL4w2s1e1.FACTORY.numberWindows());
  }}

