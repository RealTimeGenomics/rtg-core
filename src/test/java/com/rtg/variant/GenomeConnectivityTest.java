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
    final GenomeRelationships gr = new GenomeRelationships();
    gr.addGenome("father");
    gr.addGenome("mother");
    gr.addGenome("son");
    gr.addGenome("cousin");
    gr.addParentChild("father", "son");
    gr.addParentChild("mother", "son");

    final ArrayList<String> genomes = new ArrayList<>();
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
