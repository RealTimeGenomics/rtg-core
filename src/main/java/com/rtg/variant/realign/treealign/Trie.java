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

package com.rtg.variant.realign.treealign;

import java.io.IOException;
import java.io.Reader;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Trie for holding counts from alignments.
 * The trie can be constructed either by calling increment multiple times to
 * add individual byte sequences or by reading from a string. The <code>toString()</code>
 * on Trie can be used to generate a string readable this way.
 */
public class Trie extends IntegralAbstract {

  /**
   * @param in reader containing textual representation of the trie.
   * @return a trie constructed from in.
   * @throws IOException whenever.
   */
  public static Trie readTrie(final Reader in) throws IOException {
    return new ReadTrie(in).readTrie();
  }

  private static final String INDENT = "   ";

  /** A dummy array of children of length 1. Surprisingly this simplifies the code significantly. */
  final SubTrie[] mRoot = new SubTrie[1];

  private boolean mFrozen = false;

  /**
   * Increment the counts in the trie using the bytes values from start to end.
   * @param bytes values of nucleotides on range 0=N 1=A ... 4=T
   * @param start first position used in bytes (0 based)
   * @param end last position in bytes (0 based, exclusive).
   */
  public void increment(final byte[] bytes, final int start, final int end) {
    assert !mFrozen;
    increment(mRoot, 0, bytes, start, end);
  }

  /**
   * Compute probabilities from counts and freeze the counts.
   * @param priorStop prior probability that will stop at any particular position.
   */
  void freeze(final double priorStop) {
    assert !mFrozen;
    mFrozen = true;
    final SubTrie root = mRoot[0];
    if (root != null) {
      root.freeze(priorStop, 1.0);
    }
  }

  /**
   * @return true iff the counts have been frozen.
   */
  boolean frozen() {
    return mFrozen;
  }

  static void increment(final SubTrie[] parent, final int child, final byte[] bytes, final int start, final int end) {
    final SubTrie node;
    if (parent[child] == null) {
      node = new SubTrie();
      parent[child] = node;
    } else {
      node = parent[child];
    }
    node.increment(bytes, start, end);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final double total = globalIntegrity(mRoot[0]);
    if (mFrozen) {
      Exam.assertEquals(1.0, total, 1e-7);
    }
    return true;
  }

  private double globalIntegrity(final SubTrie node) {
    if (node == null) {
      return 0.0;
    }
    node.globalIntegrity();
    double sum = node.mStopProbability;
    for (int i = 0; i < node.mChildren.length; ++i) {
      sum += globalIntegrity(node.mChildren[i]);
    }
    return sum;
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(1, mRoot.length);
    return true;
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (this == obj) {
      return true;
    }
    final Trie that = (Trie) obj;
    final SubTrie thisCh = this.mRoot[0];
    final SubTrie thatCh = that.mRoot[0];
    if (thisCh == null && thatCh == null) {
      return true;
    }
    if (!(thisCh != null && thatCh != null)) {
      return false;
    }
    return thisCh.equals(thatCh);
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public void toString(final StringBuilder sb) {
    toString(sb, mRoot[0], 0);
  }

  private void toString(final StringBuilder sb, final SubTrie node, final int indent) {
    for (int i = 0; i < indent; ++i) {
      sb.append(INDENT);
    }
    if (node == null) {
      sb.append("0 0");
      sb.append(StringUtils.LS);
      return;
    }
    sb.append(node.mTotalCount).append(" ").append(node.mStopCount);
    if (mFrozen) {
      sb.append(" ").append(Utils.realFormat(node.mStopProbability, 3));
    }
    sb.append(StringUtils.LS);
    final int cont = node.mTotalCount - node.mStopCount;
    assert cont >= 0;
    if (cont > 0) {
      for (int i = 0; i < node.mChildren.length; ++i) {
        toString(sb, node.mChildren[i], indent + 1);
      }
    }
  }


}
