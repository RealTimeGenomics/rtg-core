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

package com.rtg.blacklist;

import java.io.IOException;
import java.util.HashMap;
import java.util.Random;

import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

public class HashCounterTest extends TestCase {

  @Override
  protected void setUp() throws Exception {
    super.setUp();
    Diagnostic.setLogStream();
  }

  public void test() throws IOException {
    check(2000, 18, 10);
    check(500, 24, 10);
    check(10000, 44, 10);
  }

  private void check(int length, int keyBits, int maxRepeat) throws IOException {
    final HashCounter hc = new HashCounter(length, keyBits, maxRepeat);
    final long seed = System.currentTimeMillis();
    final Random r = new Random(seed);

    //have set of roughly 80% of length is unique hashes
    final long[] keys = new long[(int) (length * 100L / 80)];
    for (int i = 0; i < keys.length; ++i) {
      keys[i] = BinaryMatrixTest.random(r, keyBits);
    }

    final long[] data = new long[length];
    //randomly fill
    for (int i = 0; i < length; ++i) {
      data[i] = keys[r.nextInt(keys.length)];
    }

    final int numThreads = 4;
    final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "HashCounterTest", false);
    try {
      for (int i = 0; i < numThreads; ++i) {
        final int threadNum = i;
        stp.execute(new IORunnable() {
          @Override
          public void run() throws IOException {
            for (int j = threadNum; j < data.length; j += numThreads) {
              try {
                hc.increment(data[j]);
              } catch (HashCounter.TooManyCollisionsException e) {
                throw new RuntimeException(e);
              }
            }
          }
        });
      }
    } finally {
      stp.terminate();
    }

    final HashMap<Long, Integer> expected = new HashMap<>();
    for (final long aData : data) {
      Integer current = expected.get(aData);
      if (current == null) {
        current = 0;
      }
      expected.put(aData, current + 1);
    }
    while (hc.next()) {
      final Integer exp = expected.remove(hc.getKey());
      assertNotNull("Seed=" + seed, exp);
      final int expInt = exp;
      assertEquals("Seed=" + seed, (long) expInt, hc.getCount());
    }

    assertEquals("Seed=" + seed, 0, expected.size());
  }
}