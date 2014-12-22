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

package com.rtg.scheduler.enumtime;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * @param <E> enumeration used as secondary ordering key.
 */
public class EnumTimeId<E extends Enum<E>> extends IntegralAbstract {

  private final int mChunk;

  private final E mType;

  /**
   * @param chunk identifier for the chunk (0 based).
   * @param type of the job.
   */
  public EnumTimeId(final int chunk, final E type) {
    mChunk = chunk;
    mType = type;
  }

  /**
   * @return the time stamp.
   */
  public int time() {
    return mChunk;
  }

  /**
   * @return the type (an enumeration that is secondary to the time stamp in the ordering).
   */
  public  E type() {
    return mType;
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null || !(obj instanceof EnumTimeId)) {
      return false;
    }

    final EnumTimeId<?> e = (EnumTimeId<?>) obj;
    return e.mChunk == this.mChunk && e.mType.equals(this.mType);
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mType.ordinal(), mChunk);
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mChunk).append(":").append(mType);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mType);
    return true;
  }
}
