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
package com.rtg.sam;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Vector;

import com.rtg.util.IORunnable;
import com.rtg.util.Populator;
import com.rtg.util.PopulatorFactory;
import com.rtg.util.ProgramState;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.io.Partition;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * A multiple threaded <code>SamFileIterator</code> which uses a <code>SortedSet</code> of
 * <code>SamMultifileIteratorRunner</code> which are run in a <code>SimpleThreadPool</code>
 * to get the parsing and validation done in the threads and then
 * allow the object using this iterator to get the current first
 * <code>SAMRecord</code>.
 *
 * TODO testing has shown that the current implementation is not
 * IO bound, context-switch bound or CPU bound so there are still
 * some possible improvements to be found.
 * @param <T> record type
 * Code review done: 28 October 2010.
 */
public final class ThreadedMultifileIterator<T> implements RecordIterator<T> {

  private static final int DEFAULT_PACKET_SIZE = 100;
  private boolean mIsClosed = false;

  private SimpleThreadPool mPool;
  /** Every member has a valid current record and its position is after the current position.  WARNING: never mutate a member of the set.  Instead, remove and re-add. */
  private final PriorityQueue<MultifileIteratorRunner<T>> mLeftmostRecordIteratorSet = new PriorityQueue<>();

  /** All the runners.  Queue used for sharpen reasons */
  private final Queue<MultifileIteratorRunner<T>> mOriginals = new LinkedList<>();
  /** All headers are identical, except for read group information. */
  private final SAMFileHeader mHeader;

  /**
   * Constructor for people wanting to run single-threaded without filtering or CRAM.
   *
   * @param files SAM/BAM files
   * @param populatorFactory populator factory
   * @throws IOException if an IO error occurs
   */
  public ThreadedMultifileIterator(Collection<File> files, PopulatorFactory<T> populatorFactory) throws IOException {
    this(new SamReadingContext(files, 1, new SamFilterParams.SamFilterParamsBuilder().create(), SamUtils.getUberHeader(files), null), populatorFactory);
  }

  /**
   * Constructor
   *
   * @param context the SAM reading context
   * @param populatorFactory populator factory
   * @throws IOException if an IO error occurs
   */
  public ThreadedMultifileIterator(final SamReadingContext context, final PopulatorFactory<T> populatorFactory) throws IOException {
    if (context.files().isEmpty()) {
      throw new IllegalArgumentException("File list is empty");
    }
    if (context.numThreads() <= 0) {
      throw new IllegalArgumentException("Illegal number of threads: " + context.numThreads());
    }

    final Collection<File> nonEmptyFiles = new ArrayList<>();
    for (File f : context.files()) {
      if (!f.isFile()) {
        nonEmptyFiles.add(f); // Pipes etc.
      } else if (f.exists() && f.length() > 0) {
        nonEmptyFiles.add(f);
      }
    }

    final int actualNumThreads = Math.min(nonEmptyFiles.size(), context.numThreads());

    final List<List<File>> fileLists = Partition.partition(actualNumThreads, nonEmptyFiles);
    mPool = new SimpleThreadPool(actualNumThreads, "ThreadedMultifileIterator", true);
    int id = 0;
    try {
      // Construction of MultifileIteratorRunner can be time consuming if doing a lot of multi-region offset querying.
      // So construct these in parallel
      final SimpleThreadPool stp = new SimpleThreadPool(actualNumThreads, "ThreadedMultifileIterator-init", false);
      final List<MultifileIteratorRunner<T>> runners = new Vector<>(); // Vector since we need thread-safe addition
      for (final List<File> subFiles : fileLists) {
        final int runnerId = id++;
        final SamReadingContext subcontext = new SamReadingContext(subFiles, 1, context.filterParams(), context.header(), context.reference(), context.referenceRanges());
        final Populator<T> populator = populatorFactory.populator();

        stp.execute(new IORunnable() {
          @Override
          public void run() throws IOException {
            runners.add(new MultifileIteratorRunner<>(subcontext, populator, runnerId, DEFAULT_PACKET_SIZE));
          }
        });
      }
      stp.terminate();
      for (final MultifileIteratorRunner<T> runner : runners) {
        mPool.execute(runner);
        mOriginals.add(runner); // So we close it
      }

      for (final MultifileIteratorRunner<T> runner : mOriginals) {
        if (runner.hasNext()) {
          mLeftmostRecordIteratorSet.add(runner);
        }
      }
    } catch (final RuntimeException | IOException e) {
      try {
        close();
      } catch (Throwable t) {
        t.addSuppressed(e);
        throw t;
      }
      throw e;
    }

    mHeader = context.header();
  }


  @Override
  public SAMFileHeader header() {
    return mHeader;
  }

  @Override
  public boolean hasNext() {
    return !mLeftmostRecordIteratorSet.isEmpty();
  }

  @Override
  public T next() {
    final MultifileIteratorRunner<T> first = mLeftmostRecordIteratorSet.poll();
    try {
      final T next = first.next();
      if (first.hasNext()) {
        mLeftmostRecordIteratorSet.add(first);
      }
      return next;
    } catch (ProgramState.SlimAbortException e) {
      try {
        close(); // This should throw the exception which was the cause of the abort
      } catch (IOException e2) { // Which, if it is an IOException, we cannot throw from Iterator, so wrap in RuntimeIOException, sigh.
        throw new RuntimeIOException(e2.getMessage(), e2);
      }
      throw e; // Fallback
    }
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }

  @Override
  public long getTotalNucleotides() {
    long count = 0;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      count += smfir.getTotalNucleotides();
    }
    return count;
  }

  /**
   * sums over all sub-iterators to get the count
   * Will only function if this iterator is closed
   * @return the sum of all invalid counts
   */
  @Override
  public long getInvalidRecordsCount() {
    long count = 0;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      count += smfir.getInvalidRecordsCount();
    }
    return count;
  }

  @Override
  public long getOutputRecordsCount() {
    long count = 0;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      count += smfir.getOutputRecordsCount();
    }
    return count;
  }

  @Override
  public long getDuplicateRecordsCount() {
    long count = 0;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      count += smfir.getDuplicateRecordsCount();
    }
    return count;
  }

  @Override
  public long getFilteredRecordsCount() {
    long count = 0;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      count += smfir.getFilteredRecordsCount();
    }
    return count;
  }

  @Override
  public long getTotalRecordsCount() {
    long count = 0;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      count += smfir.getTotalRecordsCount();
    }
    return count;
  }

  @Override
  public void close() throws IOException {
    if (mIsClosed) {
      return;
    }
    IOException problem = null;
    for (final MultifileIteratorRunner<T> smfir : mOriginals) {
      try {
        smfir.close();
      } catch (final IOException e) {
        problem = e;
      }
    }
    if (mPool != null) {
      final SimpleThreadPool s = mPool;
      mPool = null;
      s.terminate();
    }
    mIsClosed = true;
    if (problem != null) {
      throw problem;
    }
  }

}
