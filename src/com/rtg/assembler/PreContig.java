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

package com.rtg.assembler;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Contig;
import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;
import com.rtg.variant.dna.DNARange;

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
