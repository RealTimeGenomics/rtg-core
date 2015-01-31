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
import java.util.Collection;
import java.util.LinkedList;
import java.util.Queue;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.launcher.GlobalFlags;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 */
class MultifileIterator implements RecordIterator<SAMRecord> {

  static final boolean FALLBACK = GlobalFlags.isSet(GlobalFlags.SAM_ALLOW_FALLBACK_FOR_NON_INDEXED_REGIONS);

  private static final double MEGABYTE = 1 << 20;

  private final SortedSet<SamFileAndRecord> mLeftmostPriorityQueue = new TreeSet<>();
  private final Queue<SamFileAndRecord> mOriginals = new LinkedList<>();
  private final SAMFileHeader mHeader;
  private final SamFilterParams mFilterParams;

  private final long mTotalRawInputLength;
  private long mTime = 0;

  private boolean mIsClosed = false;
  private SAMRecord mNextRecordToReturn = null;


  /**
   * Constructor (deprecated)
   *
   * @param files SAM/BAM files
   * @param filterParams filter parameters. may be null for no filtering
   * @param headerOverride use this header instead of one present in file.
   * @throws IOException if an IO error occurs
   */
  public MultifileIterator(final Collection<File> files, SamFilterParams filterParams, SAMFileHeader headerOverride) throws IOException {
    this(new SamReadingContext(files, 1, filterParams, headerOverride));
  }

  /**
   * Constructor
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
    int fileCount = 0;
    long totalFileLength = 0;
    SamFileAndRecord first = null;
    try {
      for (final File file : context.files()) {
        totalFileLength += file.length();
        final RecordIterator<SAMRecord> adaptor;    // Initial source of (possibly region-restricted) SAMRecords
        final boolean streamOk = context.referenceRanges() == null || context.referenceRanges().allAvailable();
        if (file.isFile() && (!FALLBACK || streamOk || isIndexed(file))) {
          adaptor = new SamClosedFileReader(file, context.referenceRanges(), context.header());
        } else { // Fall back to SamFileAndRecord for non-file (i.e. pipes)
          Diagnostic.userLog("Using fallback for non-file or non-indexed SAM source: " + file.toString());
          adaptor = new SamFileReaderAdaptor(new SAMFileReader(FileUtils.createInputStream(file, true)), context.referenceRanges());
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

  static boolean isIndexed(File file) throws IOException {
    if (SamUtils.isBAMFile(file)) {
      return BamIndexer.indexFileName(file).exists() || BamIndexer.secondaryIndexFileName(file).exists();
    } else {
      return TabixIndexer.indexFileName(file).exists();
    }
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
  private long mInvalidRecords = 0;  // Removed here as invalid according to criteria of the calling program
  private long mOutputRecords = 0; // Made available to next level
  private static volatile long sReportedInvalidRecords = 0;

  // XXX More things than just the variant caller use this class and they might not want these being dropped
  // we should have some way to disable this notion of what is valid (or only enable it when doing variant calling)
  private boolean validRecord(final SAMRecord rec) {
    final Integer nh = SamUtils.getNHOrIH(rec);
    return (nh == null || nh > 0) && (rec.getReadUnmappedFlag() || rec.getAlignmentStart() > 0);
  }

  private boolean filterRecord(final SAMRecord r) {
    if (mFilterParams == null) {
      mOutputRecords++;
      return false;
    }
    if (!validRecord(r)) {
      mInvalidRecords++;
      if (++sReportedInvalidRecords <= 5) {
        Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING1, r.getSAMString().trim());
        Diagnostic.userLog("Invalid record: " + r.getSAMString().trim());
      }
      return true;
    } else if (!DefaultSamFilter.acceptRecord(mFilterParams, r)) {
      mFilteredRecords++;
      return true;
    }
    mOutputRecords++;
    return false;
  }

  @Override
  public long getInvalidRecordsCount() {
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getInvalidRecordsCount();
    }
    return count + mInvalidRecords;
  }

  @Override
  public long getFilteredRecordsCount() {
    long count = 0;
    for (final SamFileAndRecord sfr : mOriginals) {
      count += sfr.getInvalidRecordsCount();
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
      final SamFileAndRecord first = mLeftmostPriorityQueue.first();
      // Need to remove and re-add to priority queue to ensure correct queue behaviour
      mLeftmostPriorityQueue.remove(first);
      final SAMRecord next = first.next();
      assert next != null;
      if (first.hasNext()) {
        mLeftmostPriorityQueue.add(first);
      }
      if (!filterRecord(next)) {
        mNextRecordToReturn = next;
        hasNext = true;
        break;
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
