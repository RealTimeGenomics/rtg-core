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

package com.rtg.assembler;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Contig;
import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;
import com.rtg.mode.DNARange;

/**
*/
@TestClass({"com.rtg.assembler.PreContigTest", "com.rtg.assembler.HashMapDeBruijnGraphTest"})
class PreContig implements Contig {
  StringBuilder mContig = new StringBuilder();
  final long mId;
  int mKmerCount;

  PreContig(long id, Kmer k, int count) {
    for (int i = 0; i < k.length(); ++i) {
      mContig.append(DNA.valueChars()[k.nt(i)]);
    }
    mId = id;
    mKmerCount = count;
  }
  void extend(boolean direction, Kmer hash, int count) {
    if (direction) {
      mContig.append(DNA.valueChars()[hash.nt(hash.length() - 1)]);
    } else {
      mContig.insert(0, DNA.valueChars()[hash.nt(0)]);
    }
    mKmerCount += count;
  }
  @Override
  public int length() {
    return mContig.length();
  }

  @Override
  public byte nt(int index) {
    return DNARange.RANGE.valueOf(mContig.charAt(index));
  }

  @Override
  public String toString() {
    return "PreContig: id=" + mId + " sequence:" + mContig.toString() + StringUtils.LS;
  }

}
