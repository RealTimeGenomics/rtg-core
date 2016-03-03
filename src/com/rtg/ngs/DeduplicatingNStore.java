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
package com.rtg.ngs;

import com.rtg.util.BitPack2IntoLong;
import com.rtg.util.array.intindex.IntCreate;
import com.rtg.util.array.intindex.IntIndex;
import com.rtg.util.array.longindex.LongCreate;
import com.rtg.util.array.longindex.LongIndex;

/**
 * Keeps up to <code>n</code> records per read. If <code>n</code> is exceeded then all records for the read are dropped. Removes duplicates
 * automatically, so they don't count towards <code>n</code>.
 */
public class DeduplicatingNStore implements UptoNStore {
  private static final int FRAME_FIELD_ID = 0;
  private static final int TEMPLATE_NAME_AND_POSITION_FIELD_ID = 1;
  private static final int BITS_FOR_FRAME = 2;
  private static final int BITS_FOR_TEMPLATE_NAME_AND_POSITION_FIELD_ID = 64 - BITS_FOR_FRAME;
  final BitPack2IntoLong mBitPackHelperRecord = new BitPack2IntoLong(BITS_FOR_FRAME, BITS_FOR_TEMPLATE_NAME_AND_POSITION_FIELD_ID);

  private final LongIndex mTopNRes;
  private final IntIndex mResultCounts;
  private final int mN;
  private final long mNumTemplateSeqs;
  private final long mPositionOffset;

  /**
   * @param numReads number of reads (* 2 for paired end)
   * @param numTemplateSeqs number of template sequences
   * @param n number of results to keep per read
   * @param templateMaxLength size of largest template sequence
   * @param readMaxLength maximum read length, used to determine how far off template we can go
   */
  public DeduplicatingNStore(final long numReads, final long numTemplateSeqs, final int n, final long templateMaxLength, long readMaxLength) {
    final long length = numReads * n;

    mPositionOffset = readMaxLength;
    final long product = numTemplateSeqs * (templateMaxLength + mPositionOffset * 2);
    // Try and be a little bit careful about detecting overflow ..., also reserve 1 position
    if (product + 1 >= 1L << BITS_FOR_TEMPLATE_NAME_AND_POSITION_FIELD_ID || numTemplateSeqs == 0 || product / numTemplateSeqs != templateMaxLength + mPositionOffset * 2) {
      throw new RuntimeException("Number of sequences * maximum sequence length exceeds implementation limits.");
    }
    mN = n;
    mTopNRes = LongCreate.createIndex(length);
    mResultCounts = IntCreate.createIndex(numReads);
    mNumTemplateSeqs = numTemplateSeqs;
  }

  @Override
  public void process(long templateId, boolean reverse, int encodedReadId, int tStart, int scoreIndel) {
    final long resultIndex = (long) encodedReadId * mN;
    final int currentCount = mResultCounts.getInt(encodedReadId);
    if (currentCount > mN) {
      return;
    }
    final long storePos = (long) tStart + mPositionOffset;
    final long packed = mBitPackHelperRecord.packValues(reverse ? 1 : 0, storePos * mNumTemplateSeqs + templateId);
    for (long index = resultIndex; index < resultIndex + currentCount; index++) {
      final long val = mTopNRes.get(index);
      if (val == packed) {
        return; //duplicate
      } else if (val > packed) {
        //insert
        mResultCounts.set(encodedReadId, currentCount + 1);
        if (currentCount == mN) {
          return;
        }
        for (long i = resultIndex + currentCount - 1; i >= index; i--) {
          mTopNRes.set(i + 1, mTopNRes.get(i));
        }
        mTopNRes.set(index, packed);
        return;
      }
    }
    mResultCounts.set(encodedReadId, currentCount + 1);
    if (currentCount == mN) {
      return;
    }
    mTopNRes.set(resultIndex + currentCount, packed);
  }

  @Override
  public void setResults(MatchResult results, int encodedReadId) {
    final long resultIndex = (long) encodedReadId * mN;
    final long count = mResultCounts.get(encodedReadId);
    if (count > mN) {
      return;
    }
    for (int i = 0; i < count; i++) {
      final long result = mTopNRes.get(resultIndex + i);
      final long templateAndPos = mBitPackHelperRecord.getField(TEMPLATE_NAME_AND_POSITION_FIELD_ID, result);
      final int templateId = (int) (templateAndPos % mNumTemplateSeqs);
      final long storePos = templateAndPos / mNumTemplateSeqs;
      final int position = (int) (storePos - mPositionOffset);
      results.addMatchResult(templateId, position, encodedReadId, mBitPackHelperRecord.getField(FRAME_FIELD_ID, result) == 1L);
    }
  }

  @Override
  public String histogram() {
    return "";
  }

  @Override
  public String toString() {
    return "DeduplicatingNStore";
  }
}
