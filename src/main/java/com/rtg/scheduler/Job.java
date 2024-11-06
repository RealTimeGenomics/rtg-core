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

package com.rtg.scheduler;

import java.io.IOException;

import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

import htsjdk.samtools.SAMException;

/**
 * Has a unique id and can execute some code.
 * @param <J> the type of the job identifiers.
 */
public abstract class Job<J extends JobId<? extends JobId<?>>> extends IntegralAbstract {
  private final J mId;

  /**
   * @param id identifier.
   */
  public Job(J id) {
    mId = id;
  }

  /**
   * @return the unique identifier for this job.
   */
  public J id() {
    return mId;
  }

  /**
   * Execute the work for this job and return the result.
   * Intercepts exceptions and says what the job is that is being run.
   * @return the result (may be null).
   * @throws IOException whenever.
   */
  public final Result runCatch() throws IOException {
    try {
      return run();
    } catch (final NoTalkbackSlimException | SAMException e) {
      throw e;
    } catch (final RuntimeException e) {
      throw new RuntimeException("Job " + toString(), e);
    }
  }

  /**
   * Execute the work for this job and return the result.
   * @return the result (may be null).
   * @throws IOException whenever.
   */
  protected abstract Result run() throws IOException;

  @Override
  public String toString() {
    return id().toString();
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mId);
    return true;
  }
}
