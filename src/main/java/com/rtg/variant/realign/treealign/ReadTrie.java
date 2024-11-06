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
