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
 * An ordered set that has just  the functionality used by the scheduler.
 * @param <T> class of items in queue.
 */
public interface EventList<T> {

  /**
   * @param lookAhead the look ahead object
   * @return the first item in the queue.
   */
  T next(LookAhead lookAhead);

  /**
   * Add to the queue and insert in appropriate place.
   * @param arg item to be added.
   */
  void add(T arg);

  /**
   * @param arg item to be checked.
   * @return true iff id is contained in the queue.
   */
  boolean contains(T arg);
}
