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
package com.rtg.index;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.index.params.CreateParams;
import com.rtg.ngs.NgsParams;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * This class encapsulates a collection of Indexes and provides a way to create and freeze them.
 */
public class IndexSet {
  private final Index[] mIndexes;

  /**
   * Create an index set that encapsulates the provided indexes
   * @param indexes an array of pre-created indexes
   */
  public IndexSet(Index[] indexes) {
    mIndexes = indexes;
  }

  /**
   * Construct a collection of indexes as specified by the params objects.
   * This will create the indexes in a multi-threaded fashion.
   * @param params mapping parameters
   * @param indexParams relevant index creation params
   * @param windows the number of windows used by your hash function
   * @throws IOException if the multi-threading doesn't work out
   */
  public IndexSet(final NgsParams params, final CreateParams indexParams, int windows) throws IOException {
    mIndexes = new Index[windows];
    final Integer numberThreads = params.numberThreads();
    final SimpleThreadPool pool = new SimpleThreadPool(numberThreads, "CreateIndex", true);
    pool.enableBasicProgress(mIndexes.length);
    for (int i = 0; i < mIndexes.length; ++i) {
      pool.execute(new CreateRunnable(mIndexes, i, indexParams, params));
    }
    pool.terminate();
  }

  /**
   * @return the number of indexes encapsulated by this set
   */
  public int size() {
    return mIndexes.length;
  }

  /**
   * @param i which index would you like?
   * @return the index in position i of the set
   */
  public Index get(int i) {
    return mIndexes[i];
  }

  /**
   * Performs a multi-threaded freeze operation on all the indexes in the set.
   * @param numberThreads how many threads to use.
   * @throws IOException should the multi-threading fall over.
   */
  public void freeze(int numberThreads) throws IOException {
    final SimpleThreadPool pool = new SimpleThreadPool(numberThreads, "BuildFreeze", true);
    pool.enableBasicProgress(mIndexes.length);
    for (int i = 0; i < mIndexes.length; ++i) {
      pool.execute(new FreezeRunnable(mIndexes[i], i));
    }
    pool.terminate();
  }

  private static class CreateRunnable implements IORunnable {
    private final Index[] mIndexes;
    private final int mId;
    private final CreateParams mIndexParams;
    private final NgsParams mParams;

    CreateRunnable(final Index[] indexes, final int id, final CreateParams indexParams, final NgsParams params) {
      mIndexes = indexes;
      mId = id;
      mIndexParams = indexParams;
      mParams = params;
    }

    @Override
    public void run() {
      Diagnostic.userLog("Start create job " + mId);
      // If index hit caching is enabled this is created when threadClone is invoked rather than now.
      mIndexes[mId] = IndexUtils.createIndex(mIndexParams, mParams.indexFilter().threadClone(), mParams.numberThreads());
      Diagnostic.userLog("Finish create job " + mId);
    }
  }

  private static class FreezeRunnable implements IORunnable {
    private final Index mIndex;
    private final int mId;

    FreezeRunnable(final Index index, final int id) {
      mIndex = index;
      mId = id;
    }

    @Override
    public void run() {
      //build the internal indexes which allow searching to be done
      Diagnostic.userLog("Start freeze job " + mId);
      mIndex.freeze();
      Diagnostic.userLog("Finish freeze job " + mId);
      //System.err.println("index[" + i + "]" + LS + indexes[i]);
      Diagnostic.userLog("Index[" + mId + "] statistics " + LS + mIndex.infoString());
    }
  }
}
