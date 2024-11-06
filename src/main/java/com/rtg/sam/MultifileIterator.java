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
package com.rtg.sam;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;

/**
 */
public class MultifileIterator implements RecordIterator<SAMRecord> {

  static final boolean FALLBACK = GlobalFlags.isSet(CoreGlobalFlags.SAM_ALLOW_FALLBACK_FOR_NON_INDEXED_REGIONS);

  private static final double MEGABYTE = 1 << 20;

  private final PriorityQueue<SamFileAndRecord> mLeftmostPriorityQueue = new PriorityQueue<>();
  private final Queue<SamFileAndRecord> mOriginals = new LinkedList<>();
  private final SAMFileHeader mHeader;
  private final SamFilter mFilter;

  private final long mTotalRawInputLength;
  private long mTime = 0;

  private boolean mIsClosed = false;
  private SAMRecord mNextRecordToReturn = null;
  private String mNextRecordSourceName = null;

  /**
   * Constructor.
   *
   * @param context the SAM reading context
   * @throws IOException if an IO error occurs
   */
  public MultifileIterator(SamReadingContext context) throws IOException {
    if (context.header() == null) {
      throw new NullPointerException();
    }
    if (context.files().isEmpty()) {
      throw new IllegalArgumentException("File list is empty!");
    }
    final SamFilterParams filterParams = context.filterParams();
    mFilter = filterParams == null ? new NoneFilter() : new DefaultSamFilter(filterParams);
    int fileCount = 0;
    long totalFileLength = 0;
    SamFileAndRecord first = null;
    try {
      for (final File file : context.files()) {
        totalFileLength += file.length();
        final RecordIterator<SAMRecord> adaptor;    // Initial source of (possibly region-restricted) SAMRecords
        try {
          final boolean streamOk = context.referenceRanges() == null || context.referenceRanges().allAvailable();
          if (file.isFile() && (!FALLBACK || streamOk || SamUtils.isIndexed(file))) {
            adaptor = new SamClosedFileReader(file, context.referenceRanges(), context.reference(), context.header());
          } else { // Fall back to SamFileAndRecord for non-file (i.e. pipes)
            Diagnostic.userLog("Using fallback for non-file or non-indexed SAM source: " + file);
            adaptor = new SamFileReaderAdaptor(SamUtils.makeSamReader(FileUtils.createInputStream(file, true), context.reference(), context.header()), context.referenceRanges());
          }
        } catch (IOException | RuntimeIOException | RuntimeEOFException e) {
          throw new IOException(file.toString() + ": " + e.getMessage(), e);
        }
        final SamFileAndRecord sfr = new SamFileAndRecord(file.getPath(), fileCount++, adaptor); // Adds invalid record skipping and input source tracking
        if (first == null) {
          first = sfr;
        } else if (!SamUtils.checkHeaderDictionary(first.header(), sfr.header())) {
          mOriginals.add(sfr); // So we close it
          Diagnostic.warning(WarningType.SAM_INCOMPATIBLE_HEADERS, first.name(), sfr.name());
          throw new NoTalkbackSlimException(ErrorType.SAM_INCOMPATIBLE_HEADER_ERROR, "1");
        }

        mOriginals.add(sfr);
        if (sfr.hasNext()) {
          mLeftmostPriorityQueue.add(sfr);
        }
      }
    } catch (final RuntimeException e) {
      close();
      throw e;
    }

    mTotalRawInputLength = totalFileLength; // so we can report MB/s at the end
    mHeader = context.header();
  }

  @Override
  public final void close() throws IOException {
    if (!mIsClosed) {
      mIsClosed = true;
      //System.err.println("numFiles " + mOriginals.size());
      for (final SamFileAndRecord sfr : mOriginals) {
        sfr.close();
      }
      final double seconds = mTime / 1000000000.0;
      Diagnostic.developerLog("SamMultifileIterator IO time = " + Utils.realFormat(seconds, 3) + " s, rate = " + Utils.realFormat(mTotalRawInputLength / MEGABYTE / seconds, 3) + " MB/s");
    }
  }

  private long mFilteredRecords = 0; // Removed here due to user criteria
  private long mOutputRecords = 0; // Made available to next level

  @Override
  public long getInvalidRecordsCount() {
    return mOriginals.stream().mapToLong(RecordCounter::getInvalidRecordsCount).sum();
  }

  @Override
  public long getFilteredRecordsCount() {
    return mOriginals.stream().mapToLong(RecordCounter::getFilteredRecordsCount).sum() + mFilteredRecords;
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mOriginals.stream().mapToLong(RecordCounter::getDuplicateRecordsCount).sum();
  }

  @Override
  public long getOverCoverageRecordsCount() {
    return mOriginals.stream().mapToLong(RecordCounter::getOverCoverageRecordsCount).sum();
  }

  @Override
  public long getTotalRecordsCount() {
    return mOriginals.stream().mapToLong(RecordCounter::getTotalRecordsCount).sum();
  }

  @Override
  public long getOutputRecordsCount() {
    return mOutputRecords;
  }

  @Override
  public long getTotalNucleotides() {
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getTotalNucleotides();
    }
    return count;
  }


  @Override
  public boolean hasNext() {
    if (mNextRecordToReturn != null) {
      return true;
    }
    final long start = System.nanoTime();
    boolean hasNext = false;
    while (!mLeftmostPriorityQueue.isEmpty()) {
      // Need to remove and re-add to priority queue to ensure correct queue behaviour
      final SamFileAndRecord first = mLeftmostPriorityQueue.poll();
      final SAMRecord next = first.next();
      assert next != null;
      if (first.hasNext()) {
        mLeftmostPriorityQueue.add(first);
      }
      if (mFilter.acceptRecord(next)) {
        ++mOutputRecords;
        mNextRecordToReturn = next;
        mNextRecordSourceName = first.name();
        hasNext = true;
        break;
      } else {
        ++mFilteredRecords;
      }
    }
    mTime += System.nanoTime() - start;
    assert hasNext || mNextRecordToReturn == null;
    return hasNext;
  }

  @Override
  public SAMRecord next() {
    final SAMRecord rec = mNextRecordToReturn;
    mNextRecordToReturn = null;
    return rec;
  }

  /**
   * Return the name of the source for the record most recently returned by <code>next()</code>.
   * @return filename
   */
  public String getSource() {
    return mNextRecordSourceName;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }

  @Override
  public SAMFileHeader header() {
    return mHeader;
  }

}
