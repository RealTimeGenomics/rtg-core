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

package com.rtg.assembler.graph.io;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import com.rtg.util.Pair;
import com.rtg.util.store.StoreDirResource;
import com.rtg.util.store.StoreDirectory;

import junit.framework.TestCase;

/**
 */
public class StoreDirMutatorTest extends TestCase {

  public void test() throws IOException {
    final StoreDirectory dir = new StoreDirResource("com/rtg/assembler/planning/resources");
    final Set<Pair<String, String>> rules = new HashSet<>();
    rules.add(new Pair<>("example", "foobar"));
    rules.add(new Pair<>("\t-1", "42"));
    final StoreDirectory mutant = new StoreDirMutator(dir, rules);
    //    for (final String child : mutant.children()) {
    //      System.err.println("=== " + child);
    //      System.err.println(mutant.child(child).content());
    //    }
    final String header = mutant.child("header.txt").content();
    assertTrue(header.contains("an foobar version"));
    final String path = mutant.child("path.1.txt").content();
    //System.err.println(path);
    assertTrue(path.contains("contigs\t+242\t"));

    assertEquals("[contig.1.fa, contig.2.fa, header.txt, manifest.txt, path.1.txt]", mutant.children().toString());
    assertTrue(mutant.childExists("header.txt"));
    assertFalse(mutant.childExists("foobar.txt"));
  }
}
