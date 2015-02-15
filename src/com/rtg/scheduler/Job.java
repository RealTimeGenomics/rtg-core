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

package com.rtg.scheduler;

import java.io.IOException;

import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

import net.sf.samtools.SAMException;

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
