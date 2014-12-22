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
      mSources = new ArrayList<>();
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
