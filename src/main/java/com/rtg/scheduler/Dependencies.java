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
