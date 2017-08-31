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
