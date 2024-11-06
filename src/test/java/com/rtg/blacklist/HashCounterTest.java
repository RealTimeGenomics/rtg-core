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
          public void run() {
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
