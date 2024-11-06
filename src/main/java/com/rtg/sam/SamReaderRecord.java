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
package com.rtg.sam;

import com.rtg.util.CompareHelper;
import com.rtg.util.MathUtils;
import com.rtg.util.Populator;

import htsjdk.samtools.SAMRecord;

/**
 * This is kind of manky. Duplicate removal comes after all files are merged, and at this stage we
 * deal with ReaderRecord implementations that have normally pulled out just the information they need for
 * their task, discarding the original <code>SAMRecord</code>. So in order to do duplicate finding on a raw <code>SAMRecord</code>,
 * we need this thin wrapper around <code>SAMRecord</code> that adds the functions expected for duplicate finding.
 */
public final class SamReaderRecord implements MateInfo, ReaderRecord<SamReaderRecord> {

  static final class SamReaderRecordPopulator implements Populator<SamReaderRecord> {
    @Override
    public SamReaderRecord overflow(int position, int length) {
      throw new UnsupportedOperationException();
    }

    @Override
    public SamReaderRecord populate(SAMRecord source) {
      return new SamReaderRecord(source);
    }
  }


  final SAMRecord mRecord;
  private SamReaderRecord mNextInChain;

  private SamReaderRecord(SAMRecord rec) {
    mRecord = rec;
  }

  @Override
  public String toString() {
    return mRecord.getSAMString();
  }

  @Override
  public SamReaderRecord chain() {
    return mNextInChain;
  }

  @Override
  public void setNextInChain(SamReaderRecord rec) {
    mNextInChain = rec;
  }

  @Override
  public int disambiguateDuplicate(SamReaderRecord rec) {
    final CompareHelper helper = new CompareHelper()
      .compare(getStart(), rec.getStart())
      .compare(inner().getCigarString(), rec.inner().getCigarString())
      .compare(compare(inner().getReadBases(), rec.inner().getReadBases()))
      .compare(compare(inner().getBaseQualities(), rec.inner().getBaseQualities()));
    return helper.result();
  }

  private static int compare(final byte[] a, final byte[] b) {
    if (a.length != b.length) {
      return a.length - b.length;
    }
    for (int k = 0; k < a.length; ++k) {
      if (a[k] != b[k]) {
        return a[k] - b[k];
      }
    }
    return 0;
  }

  private static int compare(final String a, final String b) {
    if (a != null) {
      return a.compareTo(b);
    } else if (b != null) {
      return -1;
    }
    return 0;
  }

  @Override
  public int compareTo(SamReaderRecord var) {
    final int ov = valueCompareTo(var);
    if (ov != 0) {
      return ov;
    }
    return System.identityHashCode(this) - System.identityHashCode(var);
  }

  /**
   * Compares the object just on its values. the regular <code>compareTo</code> method will only return 0 for the same instance
   * @param var the object to compare
   * @return <code>-ve</code> for this object is less than, <code>0</code> for equal
   */
  int valueCompareTo(SamReaderRecord var) {
    // call starting at leftmost sequence, leftmost position, leftmost end
    if (var == this) {
      return 0;
    }
    final int thisRef = getSequenceId();
    final int thatRef = var.getSequenceId();
    if (thisRef == -1 && thatRef != -1) {
      return 1;
    } else if (thatRef == -1 && thisRef != -1) {
      return -1;
    }
    if (thisRef < thatRef) {
      return -1;
    }
    if (thisRef > thatRef) {
      return 1;
    }

    if (var.getStart() > getStart()) {
      return -1;
    } else if (var.getStart() < getStart()) {
      return 1;
    }

    final int aLen = getLength();
    final int bLen = var.getLength();
    if (aLen > bLen) {
      return 1;
    } else if (aLen < bLen) {
      return -1;
    }

    final int c = var.inner().getCigarString().compareTo(inner().getCigarString());
    if (c != 0) {
      return c;
    }
    final int mq = var.inner().getMappingQuality() - inner().getMappingQuality();
    if (mq != 0) {
      return mq;
    }
    final int r = compare(var.inner().getReadBases(), inner().getReadBases());
    if (r != 0) {
      return r;
    }
    final int q = compare(var.inner().getBaseQualities(), inner().getBaseQualities());
    if (q != 0) {
      return q;
    }
    final int f = var.inner().getFlags() - inner().getFlags();
    if (f != 0) {
      return f;
    }
    final int vascore = MathUtils.unboxNatural(var.inner().getIntegerAttribute("AS"));
    final int ascore = MathUtils.unboxNatural(inner().getIntegerAttribute("AS"));
    final int as = vascore - ascore;
    if (as != 0) {
      return as;
    }
    final int vamb = MathUtils.unboxNatural(SamUtils.getNHOrIH(var.inner()));
    final int amb = MathUtils.unboxNatural(SamUtils.getNHOrIH(inner()));
    final int nh = vamb - amb;
    if (nh != 0) {
      return nh;
    }
    final String vss = var.inner().getStringAttribute(SamUtils.CG_SUPER_CIGAR);
    final String ss = inner().getStringAttribute(SamUtils.CG_SUPER_CIGAR);
    final int sc = compare(vss, ss);
    if (sc != 0) {
      return sc;
    }
    return 0;
  }

  @Override
  public boolean equals(final Object var) {
    return this == var;
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public int getGenome() {
    return 0;
  }

  /**
   * @return the <code>SAMRecord</code> being wrapped
   */
  public SAMRecord inner() {
    return mRecord;
  }

  @Override
  public boolean isMated() {
    return mRecord.getReadPairedFlag() && mRecord.getProperPairFlag();
  }

  @Override
  public int getFragmentLength() {
    return mRecord.getInferredInsertSize();
  }

  @Override
  public int getMateSequenceId() {
    return mRecord.getMateReferenceIndex();
  }

  @Override
  public int getSequenceId() {
    return mRecord.getReferenceIndex();
  }

  @Override
  public int getStart() {
    return mRecord.getAlignmentStart() - 1;
  }

  @Override
  public int getEnd() {
    return mRecord.getAlignmentEnd();
  }

  @Override
  public int getLength() {
    return getEnd() - getStart();
  }

}
