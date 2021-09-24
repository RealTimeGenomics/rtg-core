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

/**
 * Immutable labels for <code>Job</code>s. The methods declared are there as a reminder
 * that they really do need to be implemented and that users of <code>JobId</code> rely on them.
 * @param <T> the type used by comparable. A particular implementation type used consistently for all jobs.
 */
public interface JobId<T> extends Comparable<T> {

  /**
   * @param results to be validated.
   * @return true iff the number and types of the results are correct as arguments for a job with this identifier.
   */
  boolean validArguments(Result[] results);

  /**
   * @param result to be validated.
   * @return true iff the number and types of the objects in result are correct as the result for a job with this identifier.
   */
  boolean validResult(Result result);

  /**
   * @return the time stamp.
   */
  int time();

  @Override
  boolean equals(Object arg0);

  @Override
  int hashCode();

  @Override
  String toString();

}
