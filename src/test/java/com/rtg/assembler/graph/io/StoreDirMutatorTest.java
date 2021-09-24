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
