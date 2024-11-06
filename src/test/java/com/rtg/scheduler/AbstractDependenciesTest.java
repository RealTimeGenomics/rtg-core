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
import java.util.HashSet;
import java.util.Set;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 * @param <J> job id
 */
public abstract class AbstractDependenciesTest<J extends JobId<J>> extends TestCase {

  protected abstract Dependencies<J> getDependencies();

  public void test() {
    allIds(getDependencies(), 5);
  }

  protected LookAhead getLookAhead(int lookAhead, int delta) {
    return new LookAhead(lookAhead, delta);
  }

  protected Set<J> allIds(final Dependencies<J> dep, int lookAhead) {
    final Set<J> allIds = new HashSet<>();
    final LookAhead head = getLookAhead(lookAhead, dep.delta());
    while (true) {
      final J id = dep.next(head);
      //System.err.println("init:" + id);
      if (id == null) {
        break;
      }
      head.increment(id.time());
      allIds.add(id);
      assertTrue("" + id + ":" + dep.from(id), Util.nonNullSize(dep.from(id)) == 0);
    }
    while (true) {
      final Set<J> frontier = new HashSet<>();
      for (J next : allIds) {
        final Collection<J> to = dep.to(next);
        for (J idTo : to) {
          if (idTo != null && !allIds.contains(idTo)) {
            //System.err.println(next + ">" + idTo);
            frontier.add(idTo);
          }
        }
      }
      //System.err.println("Frontier:" + frontier);
      //System.err.println("All:" + allIds);
      if (frontier.size() == 0) {
        break;
      }
      allIds.addAll(frontier);
    }

    //check that from and to are ordered
    for (J next : allIds) {
      Exam.integrity(next);
      final Collection<J> to = dep.to(next);
      Util.checkOrder(next, to, +1);
      final Collection<J> from = dep.from(next);
      Util.checkOrder(next, from, -1);
    }
    return allIds;
  }

}
