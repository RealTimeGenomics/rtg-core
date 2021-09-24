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

package com.rtg.variant.realign.treealign;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.StringReader;

import junit.framework.TestCase;

/**
 */
public class TrieTest extends TestCase {

  //empty trie
  public void test0() throws IOException {
    final Trie trie = new Trie();
    assertFalse(trie.frozen());
    trie.globalIntegrity();
    final String str = trie.toString();
    assertEquals("0 0" + LS, str);
    assertEquals(trie, Trie.readTrie(new StringReader(str)));
  }

  //add three simple strings
  public void test() throws IOException {
    final Trie trie = new Trie();
    trie.increment(new byte[] {1, 1, 2}, 0, 3);
    trie.increment(new byte[] {1, 1, 3}, 0, 3);
    trie.increment(new byte[] {1, 1, 3}, 0, 1);
    trie.globalIntegrity();
    final String str = trie.toString();
    final String exp = ""
        + "3 0" + LS
        + "   3 1" + LS
        + "      2 0" + LS
        + "         0 0" + LS
        + "         1 1" + LS
        + "         1 1" + LS
        + "         0 0" + LS
        + "      0 0" + LS
        + "      0 0" + LS
        + "      0 0" + LS
        + "   0 0" + LS
        + "   0 0" + LS
        + "   0 0" + LS
        ;
    assertEquals(exp, str);
    assertEquals(trie, Trie.readTrie(new StringReader(str)));

    assertFalse(trie.frozen());
    trie.freeze(0.1);
    assertTrue(trie.frozen());
    trie.globalIntegrity();
    final String expProb = ""
        + "3 0 0.194" + LS
        + "   3 1 0.358" + LS
        + "      2 0 0.082" + LS
        + "         0 0" + LS
        + "         1 1 0.183" + LS
        + "         1 1 0.183" + LS
        + "         0 0" + LS
        + "      0 0" + LS
        + "      0 0" + LS
        + "      0 0" + LS
        + "   0 0" + LS
        + "   0 0" + LS
        + "   0 0" + LS
        ;
    assertEquals(expProb, trie.toString());
  }

  //add strings with edge cases - beyond end of array - 0 - 5
  public void testEdges() throws IOException {
    final Trie trie = new Trie();
    trie.increment(new byte[] {1, 1, 0}, 0, 3);
    trie.increment(new byte[] {1, 1, 5}, 0, 3);
    trie.increment(new byte[] {1, 1}, 2, 5);
    trie.globalIntegrity();
    final String str = trie.toString();
    final String exp = ""
        + "3 1" + LS
        + "   2 0" + LS
        + "      2 2" + LS
        + "      0 0" + LS
        + "      0 0" + LS
        + "      0 0" + LS
        + "   0 0" + LS
        + "   0 0" + LS
        + "   0 0" + LS
        ;
    assertEquals(exp, str);
    assertEquals(trie, Trie.readTrie(new StringReader(str)));
  }

  //special top-level cases
  public void testEqualsNull() {
    final Trie trie = new Trie();
    assertEquals(System.identityHashCode(trie), trie.hashCode());
    assertTrue(trie.equals(trie));
    assertFalse(trie.equals(null));
  }

  //different ojects but both have null root
  public void testEquals1() {
    final Trie trie0 = new Trie();
    final Trie trie1 = new Trie();
    assertTrue(trie0.equals(trie1));
  }

  //empty and non-empty
  public void testEquals2() {
    final Trie trie0 = new Trie();
    trie0.increment(new byte[] {1}, 0, 1);
    final Trie trie1 = new Trie();
    assertFalse(trie0.equals(trie1));
  }

}
