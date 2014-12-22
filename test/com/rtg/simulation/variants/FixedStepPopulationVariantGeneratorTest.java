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
package com.rtg.simulation.variants;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class FixedStepPopulationVariantGeneratorTest extends TestCase {

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
  }

  private static final String REF = ">ref" + StringUtils.LS
          + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
          + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
          + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca";

  public void testFixedStepX() throws IOException {
    final SequencesReader sr = ReaderTestUtils.getReaderDnaMemory(REF);
    FixedStepPopulationVariantGenerator fixed = new FixedStepPopulationVariantGenerator(sr, 10, new Mutator("X"), new PortableRandom(10), 0.5);
    List<PopulationVariantGenerator.PopulationVariant> variants = fixed.generatePopulation();
    assertEquals(12, variants.size());
    int e = 0;
    int i = REF.indexOf("c");
    for (; i < REF.length(); i += 10) {
      PopulationVariantGenerator.PopulationVariant var = variants.get(e);
      assertTrue(Arrays.equals(DnaUtils.encodeString(REF.substring(i, i + 1)), var.mRef));
      assertFalse(Arrays.equals(var.mAlleles[0], var.mRef));
      assertEquals(0, var.getSequenceId());
      assertEquals(e * 10, var.getStart());
      e++;
    }
  }
  public void testFixedStepHetX() throws IOException {
    final SequencesReader sr = ReaderTestUtils.getReaderDnaMemory(REF);
    FixedStepPopulationVariantGenerator fixed = new FixedStepPopulationVariantGenerator(sr, 10, new Mutator("X_Y"), new PortableRandom(118), 0.5);
    List<PopulationVariantGenerator.PopulationVariant> variants = fixed.generatePopulation();
    assertEquals(12, variants.size());
    int e = 0;
    int i = REF.indexOf("c");
    for (; i < REF.length(); i += 10) {
      PopulationVariantGenerator.PopulationVariant var = variants.get(e);
      assertTrue(Arrays.equals(DnaUtils.encodeString(REF.substring(i, i + 1)), var.mRef));
      assertTrue(!Arrays.equals(var.mAlleles[0], var.mAlleles[1]));
      assertEquals(0, var.getSequenceId());
      assertEquals(e * 10, var.getStart());
      e++;
    }
  }
  public void testFixedStepI() throws IOException {
    final SequencesReader sr = ReaderTestUtils.getReaderDnaMemory(REF);
    FixedStepPopulationVariantGenerator fixed = new FixedStepPopulationVariantGenerator(sr, 10, new Mutator("I"), new PortableRandom(10), 0.5);
    List<PopulationVariantGenerator.PopulationVariant> variants = fixed.generatePopulation();
    assertEquals(11, variants.size());
    int e = 0;
    int i = REF.indexOf("c") + 10;
    for (; i < REF.length(); i += 10) {
      PopulationVariantGenerator.PopulationVariant var = variants.get(e);
      assertTrue(Arrays.equals(DnaUtils.encodeString(REF.substring(i - 1, i)), var.mRef));
      assertFalse(Arrays.equals(var.mAlleles[0], var.mRef));
      assertEquals(0, var.getSequenceId());
      assertEquals((e + 1) * 10 - 1, var.getStart());
      e++;
    }
  }
}
