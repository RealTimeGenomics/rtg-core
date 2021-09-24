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

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

import com.rtg.assembler.graph.Graph;
import com.rtg.util.IORunnable;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.ProgramState;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public final class AsyncReadSource implements IORunnable {


  private static final int TIMEOUT = 1;
  private volatile boolean mClosed = false;
  private final ReadPairSource mSource;
  private final BlockingQueue<List<byte[]>> mQueue = new LinkedBlockingQueue<>(1024);
  final String mName;

  /**
   * Split reading of a read source out into a separate thread.
   * @param source the source to multi-thread
   * @param name A helpful name for the logs
   */
  AsyncReadSource(ReadPairSource source, String name) {
    mSource = source;
    mName = name;
  }

  public int getNumberFragments() {
    return mSource.numberFragments();
  }
  public long getNumberReads() {
    return mSource.numberReads();
  }

  /**
   * Delegate for aligner creation.
   * @param graph the graph to align against
   * @param mismatches mismatches allowed
   * @param traverse graph traversal structure
   * @return a new graph aligner
   */
  public GraphAligner aligner(Graph graph, IntegerOrPercentage mismatches, GraphTraversions traverse) {
    return mSource.aligner(graph, mismatches, traverse);
  }

  /**
   * Delegate for minimum insert size.
   * @return minimum insert size of underlying source
   */
  public int minInsertSize() {
    return mSource.minInsertSize();
  }

  /**
   * Delegate for maximum insert size.
   * @return maximum insert size of underlying source
   */
  public int maxInsertSize() {
    return mSource.maxInsertSize();
  }

  /**
   * Retrieve the underlying source.
   * @return the read source this wraps
   */
  public ReadPairSource underlying() {
    return mSource;
  }

  /**
   * Retrieve the next fragment set from the underlying set.
   * @return the fragments of the read as a list of byte arrays
   */
  public List<byte[]> nextFragments() {
    while (!mClosed) {
      try {
        final List<byte[]> fragments = mQueue.poll(TIMEOUT, TimeUnit.SECONDS);
        if (fragments != null) {
          if (fragments.isEmpty()) {
            // Put back onto queue, so that subsequent calls to nextFragments() get the sentinel
            mQueue.put(fragments);
            return null;
          } else {
            return fragments;
          }
        }
        ProgramState.checkAbort();
      } catch (InterruptedException e) {
        ProgramState.setAbort();
      }
    }
    return null;
  }

  @Override
  public void run() throws IOException {
    Diagnostic.developerLog("AsyncReadSource: " + mName + " started reading");
    List<byte[]> fragments;
    while (!mClosed && (fragments = mSource.nextFragments()) != null) {
      try {
        while (!mQueue.offer(fragments, TIMEOUT, TimeUnit.SECONDS)) {
          ProgramState.checkAbort();
        }
      } catch (InterruptedException e) {
        ProgramState.setAbort();
      }
    }
    if (!mClosed) {
      try {
        while (!mQueue.offer(Collections.emptyList(), TIMEOUT, TimeUnit.SECONDS)) {
          ProgramState.checkAbort();
        }
        Diagnostic.developerLog("AsyncReadSource: " + mName + " finished reading");
      } catch (InterruptedException e) {
        ProgramState.setAbort();
      }
    }
  }

  /**
   * Indicate that the threads shouldn't continue to use this source
   */
  public void close() {
    mClosed = true;
  }
}

