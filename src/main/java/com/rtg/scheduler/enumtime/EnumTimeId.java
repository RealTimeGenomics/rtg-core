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
