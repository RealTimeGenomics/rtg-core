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
  private final SamFilterParams mFilterParams;
  private final SamFilter mFilter;

  private final long mTotalRawInputLength;
  private long mTime = 0;

  private boolean mIsClosed = false;
  private SAMRecord mNextRecordToReturn = null;

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
    if (context.files().size() == 0) {
      throw new IllegalArgumentException("File list is empty!");
    }
    mFilterParams = context.filterParams();
    mFilter = mFilterParams == null ? new NoneFilter() : new DefaultSamFilter(mFilterParams);
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
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getInvalidRecordsCount();
    }
    return count;
  }

  @Override
  public long getFilteredRecordsCount() {
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getFilteredRecordsCount();
    }
    return count + mFilteredRecords;
  }

  @Override
  public long getDuplicateRecordsCount() {
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getDuplicateRecordsCount();
    }
    return count;
  }

  @Override
  public long getTotalRecordsCount() {
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getTotalRecordsCount();
    }
    return count;
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

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }

  @Override
  public SAMFileHeader header() {
    return mHeader;
  }

}
