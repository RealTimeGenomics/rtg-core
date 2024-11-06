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

package com.rtg.assembler;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.util.SimpleThreadPool;

/**
 * Wrap the construction of a simple thread pool for reading the reads asynchronously
 */
public class AsyncReadPool implements Closeable {
  private final SimpleThreadPool mPool;
  private final List<AsyncReadSource> mSources;

  /**
   * Create a thread pool to read the provided sources
   * @param name a name to use in the logs
   * @param sources read these asynchronously
   */
  public AsyncReadPool(String name, List<ReadPairSource> sources) {
    if (sources.size() > 0) {
      mPool = new SimpleThreadPool(sources.size(), name, true);
      mSources = new ArrayList<>(sources.size());
      int id = 0;
      for (ReadPairSource source : sources) {
        final AsyncReadSource async = new AsyncReadSource(source, name + "-" + id++);
        mSources.add(async);
        mPool.execute(async);
      }
    } else {
      mPool = null;
      mSources = Collections.emptyList();
    }
  }

  /**
   * Delegate for <code>SimpleThreadPool</code> terminate
   * @throws IOException if the underlying pool
   */
  public void terminate() throws IOException {
    if (mPool != null) {
      mPool.terminate();
    }
  }

  /**
   * Retrieve the asynchronous readers
   * @return the readers that have been spawned
   */
  public List<AsyncReadSource> sources() {
    return mSources;
  }

  @Override
  public void close() throws IOException {
    for (AsyncReadSource source : mSources) {
      source.close();
    }
    terminate();
  }
}
