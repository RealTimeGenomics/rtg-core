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

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.BitPack3IntoLong;
import com.rtg.util.StringUtils;
import com.rtg.util.array.longindex.LongCreate;
import com.rtg.util.array.longindex.LongIndex;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * <p>Implementation of separate from result writing, used in paired end writer.
 * <p>Should be modified so both <code>TopN</code> classes can use it and
 * remove details of the result so this class doesn't have to know what data
 * it is storing / ranking, just how to compare it.
 */
public class TopNImplementation implements UptoNStore {
  //Per record bit fields
  private static final int FRAME_FIELD_ID = 0;
  private static final int SCORE_FIELD_ID = 1;
  private static final int TEMPLATE_NAME_AND_POSITION_FIELD_ID = 2;

  //Per read bit fields
  private static final int RESULT_COUNT_FIELD_ID = 0;
  private static final int WORST_SCORE_FIELD_ID = 1;
  private static final int WORST_SCORE_COUNT_FIELD_ID = 2;

  private static final int BITS_FOR_RESULT_COUNT = 16;
  private static final int BITS_FOR_WORST_SCORE = 12;
  private static final int BITS_FOR_WORST_SCORE_COUNT = 64 - BITS_FOR_RESULT_COUNT - BITS_FOR_WORST_SCORE;
  private static final long MAX_SCORE_COUNT = (1L << BITS_FOR_WORST_SCORE_COUNT) - 1;
  private static final int MAX_SCORE = (1 << BITS_FOR_WORST_SCORE) - 1;

  private static final int BITS_FOR_FRAME = 1;
  private static final int BITS_FOR_SCORE = 12;
  private static final int BITS_FOR_TEMPLATE_NAME_AND_POSITION_FIELD_ID = 64 - BITS_FOR_SCORE - BITS_FOR_FRAME;
  final BitPack3IntoLong mBitPackHelperRecord = new BitPack3IntoLong(BITS_FOR_FRAME, BITS_FOR_SCORE, BITS_FOR_TEMPLATE_NAME_AND_POSITION_FIELD_ID);
  final BitPack3IntoLong mBitPackHelperRead = new BitPack3IntoLong(BITS_FOR_RESULT_COUNT, BITS_FOR_WORST_SCORE, BITS_FOR_WORST_SCORE_COUNT);

  //Reserved for a return value from insert, TopEquals constructor makes sure that not all bits can be set in result
  private static final long NO_RESULT = ~0L;

  private final LongIndex mTopNRes;
  private final LongIndex mResultCounts;

  private final int mN;
  private final long mNumTemplateSeqs;
  private final long mPositionOffset;

  /**
   * Constructor
   * @param numReads number of reads (* 2 for paired end)
   * @param numTemplateSeqs number of template sequences
   * @param n number of results to keep per read
   * @param templateMaxLength size of largest template sequence
   * @param readMaxLength maximum read length, used to determine how far off template we can go
   */
  public TopNImplementation(final long numReads, final long numTemplateSeqs, final int n, final long templateMaxLength, long readMaxLength) {
    final long length = numReads * n;
    mPositionOffset = readMaxLength;
    final long product = numTemplateSeqs * (templateMaxLength + mPositionOffset * 2);
    // Try and be a little bit careful about detecting overflow ..., also reserve 1 position
    if (product + 1 >= 1L << BITS_FOR_TEMPLATE_NAME_AND_POSITION_FIELD_ID || numTemplateSeqs == 0 || product / numTemplateSeqs != templateMaxLength + mPositionOffset * 2) {
      throw new RuntimeException("Number of sequences * maximum sequence length exceeds implementation limits.");
    }

    mTopNRes = LongCreate.createIndex(length);
    mResultCounts = LongCreate.createIndex(numReads);
    mN = n;
    mNumTemplateSeqs = numTemplateSeqs;
    Diagnostic.userLog(toString() + " statistics" + LS + infoString());
  }

  String infoString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Memory Usage\tbytes\tlength").append(LS);
    long totalBytes = 0;
    sb.append("\t\t").append(StringUtils.commas(mTopNRes.bytes())).append("\t").append(StringUtils.commas(mTopNRes.length())).append("\tHTopNRes").append(LS);
    totalBytes +=  mTopNRes.bytes();

    sb.append("\t\t").append(StringUtils.commas(mResultCounts.bytes())).append("\t").append(StringUtils.commas(mResultCounts.length())).append("\tResultCounts").append(LS);
    totalBytes +=  mResultCounts.bytes();

    sb.append("\t\t").append(StringUtils.commas(totalBytes)).append("\t\tTotal bytes").append(LS);

    return sb.toString();
  }

  /**
   * @return a printable histogram of the <code>TopN</code> data
   */
  @Override
  public String histogram() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Histogram for TopN: ").append(LS);
    final long[] hist = new long[mN + 1];
    for (int i = 0; i < mResultCounts.length(); ++i) {
      hist[(int) mBitPackHelperRead.getField(RESULT_COUNT_FIELD_ID, mResultCounts.get(i))]++;
    }
    long tot = 0;
    for (int i = 0; i < hist.length; ++i) {
      sb.append(i).append(": ").append(hist[i]).append(LS);
      if (i > 0) {
        tot += hist[i] * i;
      }
    }
    sb.append("TopN Utilisation: ").append(tot).append(" / ").append(mTopNRes.length()).append(" (").append(tot * 100L / mTopNRes.length()).append("%)");
    return sb.toString();
  }

  /**
   * Add a match
   * @param templateId id of template sequence
   * @param reverse true if on reverse frame
   * @param encodedReadId encoded read id
   * @param tStart template start position of match
   * @param scoreIndel indel score
   */
  @Override
  public void process(final long templateId, boolean reverse, final int encodedReadId, final int tStart, int scoreIndel) {

    final long readPacked = mResultCounts.get((long) encodedReadId);
    final long worstScoreCount = mBitPackHelperRead.getField(WORST_SCORE_COUNT_FIELD_ID, readPacked);
    final long worstScore = mBitPackHelperRead.getField(WORST_SCORE_FIELD_ID, readPacked);
    final long resultCount = mBitPackHelperRead.getField(RESULT_COUNT_FIELD_ID, readPacked); //currently maxs at 65535
    final int storedScore;
    if (scoreIndel > MAX_SCORE) {
      Diagnostic.warning("Score " + scoreIndel + " exceeds maximum score, truncating to " + MAX_SCORE);
      storedScore = scoreIndel & MAX_SCORE;
    } else {
      storedScore = scoreIndel;
    }

    if (worstScoreCount > 0) {
      if (storedScore > worstScore) {
        return;
      } else if (storedScore == worstScore) {
        final long newWorstScoreCount = worstScoreCount < MAX_SCORE_COUNT ? worstScoreCount + 1 : MAX_SCORE_COUNT;
        mResultCounts.set((long) encodedReadId, mBitPackHelperRead.packValues(resultCount, worstScore, newWorstScoreCount));
      }
    }

    final long resultIndex = (long) encodedReadId * mN;
    final long storePos = (long) tStart + mPositionOffset;
    final long result = mBitPackHelperRecord.packValues(reverse ? 1L : 0L, storedScore, storePos * mNumTemplateSeqs + templateId);
    final long drop = insert(storedScore, resultIndex, (int) resultCount, result);

    //read records is full and new value wasn't inserted (due to score), also have checked against worstScore already
    if (drop == NO_RESULT && resultCount == mN) {
      //straight into count array, is better than current one
      mResultCounts.set((long) encodedReadId, mBitPackHelperRead.packValues(resultCount, storedScore, 1));
      return;
    }

    final long newCountVal;
    if (resultCount == mN) {
      final int dropScore = (int) mBitPackHelperRecord.getField(SCORE_FIELD_ID, drop);

      if (dropScore < worstScore || worstScoreCount == 0) {
        newCountVal = mBitPackHelperRead.packValues(resultCount, dropScore, 1);
      } else {
        //things in top N are <= edgeScore
        final long newWorstScoreCount = worstScoreCount < MAX_SCORE_COUNT ? worstScoreCount + 1 : MAX_SCORE_COUNT;
        newCountVal = mBitPackHelperRead.packValues(resultCount, worstScore, newWorstScoreCount);
      }
    } else {
      final long newResultCount = resultCount + 1;
      assert newResultCount <= mN;
      newCountVal = mBitPackHelperRead.packValues(newResultCount, 0, 0);
    }
    mResultCounts.set((long) encodedReadId, newCountVal);
  }

  //returns a value because this is required by the original TopNOutputProcessor, which I would like
  //to eventually use a modified version of this class (once the details of results have been moved
  //out to the calling class)
  private long insert(final int insertScore, final long resultIndex, final int count, final long result) {
    //reverse insertion
    final long ret;
    long i;
    if (count == mN) {
      //check here to avoid bounds exception
      final long last = mTopNRes.get(count + resultIndex - 1);
      if (insertScore >=  mBitPackHelperRecord.getField(SCORE_FIELD_ID, last)) {
        return NO_RESULT;
      }
      i = count + resultIndex - 1;
      ret = last;
    } else {
      i = count + resultIndex;
      ret = NO_RESULT;
    }
    //find insert position, shifting as we go
    for (; i > resultIndex; --i) {
      final long current = mTopNRes.get(i - 1);
      if (insertScore >= mBitPackHelperRecord.getField(SCORE_FIELD_ID, current)) {
        break;
      }
      mTopNRes.set(i, mTopNRes.get(i - 1));
    }
    mTopNRes.set(i, result);
    return ret;
  }


  @Override
  public void setResults(MatchResult results, final int encodedReadId) {
    final long resultIndex = (long) encodedReadId * mN;
    final long packedCount = mResultCounts.get(encodedReadId);
    final long count = mBitPackHelperRead.getField(RESULT_COUNT_FIELD_ID, packedCount);
    final long worstScoreCount = mBitPackHelperRead.getField(WORST_SCORE_COUNT_FIELD_ID, packedCount);
    final long worstScore = worstScoreCount > 0 ? mBitPackHelperRead.getField(WORST_SCORE_FIELD_ID, packedCount) : Integer.MAX_VALUE;

    for (int i = 0; i < count; ++i) {
      final long result = mTopNRes.get(resultIndex + i);
      final int scoreIndel = (int) mBitPackHelperRecord.getField(SCORE_FIELD_ID, result);
      if (scoreIndel >= worstScore) {
        break;
      }
      final long templateAndPos = mBitPackHelperRecord.getField(TEMPLATE_NAME_AND_POSITION_FIELD_ID, result);
      final int templateId = (int) (templateAndPos % mNumTemplateSeqs);
      final long storePos = templateAndPos / mNumTemplateSeqs;
      final int position = (int) (storePos - mPositionOffset);
      results.addMatchResult(templateId, position, encodedReadId, mBitPackHelperRecord.getField(FRAME_FIELD_ID, result) == 1L);
    }
  }

  @Override
  public String toString() {
    return "TopNImplementation";
  }

}
