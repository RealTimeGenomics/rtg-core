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
