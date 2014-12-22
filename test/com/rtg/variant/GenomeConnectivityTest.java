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
package com.rtg.variant;

import java.util.ArrayList;

import com.rtg.relation.GenomeRelationships;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class GenomeConnectivityTest extends TestCase {

  public void test() {
    TestUtils.testEnum(GenomeConnectivity.class, "[SPARSE, DENSE]");
  }

  public void testGetConnectivity() {
    GenomeRelationships gr = new GenomeRelationships();
    gr.addGenome("father");
    gr.addGenome("mother");
    gr.addGenome("son");
    gr.addGenome("cousin");
    gr.addParentChild("father", "son");
    gr.addParentChild("mother", "son");

    ArrayList<String> genomes = new ArrayList<>();
    genomes.add("father");
    genomes.add("mother");
    genomes.add("son");

    assertEquals(GenomeConnectivity.DENSE, GenomeConnectivity.getConnectivity(genomes, gr));

    genomes.add("cousin");

    assertEquals(GenomeConnectivity.DENSE, GenomeConnectivity.getConnectivity(genomes, gr));

    gr.addGenome("cousin2");
    genomes.add("cousin2");

    assertEquals(GenomeConnectivity.SPARSE, GenomeConnectivity.getConnectivity(genomes, gr));

    genomes.remove("cousin2");

    gr.addGenome("father2");
    gr.addGenome("mother2");
    gr.addGenome("son2");
    gr.addParentChild("father2", "son2");
    gr.addParentChild("mother2", "son2");
    genomes.add("father2");
    genomes.add("mother2");
    genomes.add("son2");

    assertEquals(GenomeConnectivity.SPARSE, GenomeConnectivity.getConnectivity(genomes, gr));

    gr.addGenome("grandmother");
    gr.addParentChild("grandmother", "father2");
    gr.addParentChild("grandmother", "mother");
    genomes.add("grandmother");

    assertEquals(GenomeConnectivity.DENSE, GenomeConnectivity.getConnectivity(genomes, gr));
  }
}
