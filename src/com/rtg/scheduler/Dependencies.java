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

import java.util.Collection;

/**
 * Information about the dependencies between  <code>Job</code>s
 * using only <code>JobId</code>s.
 * This varies depending on the problem being run.
 * @param <J> the type of the job identifiers.
 */
public interface Dependencies<J extends JobId<J>> {

  /**
   * Finds all the <code>JobId</code>s that must have completed before
   * the job labeled by <code>id</code> is run.<br>
   * The job identifiers returned will be in the order of the arguments needed for the
   * job derived from <code>id</code>. (Note: implementors should be very careful to
   * obey this constraint, <code>LinkedHashSet</code> may be useful).
   * All the job identifiers will be ordered earlier than <code>id</code>.
   * @param id the unique job identifier being checked.
   * @return the <code>JobId</code>s (may not be null, may be empty).
   */
  Collection<J> from(J id);

  /**
   * Finds all the <code>JobId</code>s that are waiting on the results of executing <code>id</code>.
   * They will all be ordered after <code>id</code>.
   * @param id the unique job identifier being checked.
   * @return the <code>JobId</code>s (may not be null, may be empty).
   */
  Collection<J> to(J id);

  /**
   * Get another job which can be run.
   * Exhaustively enumerates all the possible job identifiers.
   * All these jobs should require only the <code>JobId</code> as arguments in order
   * to execute.
   * @param lookAhead the look ahead object
   * @return the job identifier. null if no more available.
   */
  J next(LookAhead lookAhead);

  /**
   * Get the maximum dependency depth.
   * @return the maximum dependency depth.
   */
  int delta();
}
