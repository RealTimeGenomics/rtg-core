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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Helper class to construct a trie from a reader.
 */
class ReadTrie extends IntegralAbstract {
  private final BufferedReader mIn;

  ReadTrie(final Reader in) {
    mIn = new BufferedReader(in);
  }

  Trie readTrie() throws IOException {
    final Trie trie = new Trie();
    read(trie.mRoot, 0);
    return trie;
  }

  void read(final SubTrie[] parent, final int child) throws IOException {
    final String line = mIn.readLine();
    if (line == null) {
      return;
    }
    final String[] split = line.split("\\s+");
    //System.err.println(line);
    //System.err.println(Arrays.toString(split));
    final int total;
    final int stop;
    if (split.length == 2) {
      total = Integer.parseInt(split[0]);
      stop = Integer.parseInt(split[1]);
    } else if (split.length == 3) {
      total = Integer.parseInt(split[1]);
      stop = Integer.parseInt(split[2]);
    } else {
      throw new RuntimeException("Incorrect number of fields:" + line);
    }
    if (total > 0) {
      final SubTrie sub = new SubTrie(total, stop);
      parent[child] = sub;
      if (total > stop) {
        for (int i = 0; i < sub.mChildren.length; ++i) {
          read(sub.mChildren, i);
        }
      }
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mIn);
    return true;
  }

}
