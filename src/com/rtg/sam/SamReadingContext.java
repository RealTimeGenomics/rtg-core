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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;

import htsjdk.samtools.SAMFileHeader;

/**
 * When we are about to read sam files, this object just stores the file list, the filter params and the uber header
 */
@TestClass("com.rtg.sam.DefaultSamFilterTest")
public class SamReadingContext {

  private final Collection<File> mFiles;
  private final SamFilterParams mParams;

  private final SAMFileHeader mHeader;
  private final int mNumThreads;

  private final ReferenceRanges mReferenceRanges;

  /**
   * Create the reading context, using information in the filter params to resolve restriction regions.
   * @param files the files to be read
   * @param numThreads the number of threads to use for reading
   * @param filterParams supplies settings involving individual record filtering and region restrictions
   * @param header the already obtained merged SAM header.
   * @throws IOException if the restriction involved reading a BED file that could not be read
   */
  public SamReadingContext(Collection<File> files, int numThreads, SamFilterParams filterParams, SAMFileHeader header) throws IOException {
    this(files, numThreads, filterParams, header,  SamRangeUtils.createReferenceRanges(header, filterParams));
  }

  /**
   * Create the reading context, using information in the filter params to resolve restriction regions.
   * @param files the files to be read
   * @param numThreads the number of threads to use for reading
   * @param filterParams supplies settings involving individual record filtering and region restrictions
   * @param header the already obtained merged SAM header.
   * @param referenceRanges the already resolved reference ranges to read over.
   */
  public SamReadingContext(Collection<File> files, int numThreads, SamFilterParams filterParams, SAMFileHeader header, ReferenceRanges referenceRanges) {
    if (header == null) {
      throw new NullPointerException();
    }
    mFiles = files;
    mParams = filterParams;
    mHeader = header;
    mNumThreads = numThreads;
    mReferenceRanges = referenceRanges;
  }

  /**
   * @return the SAM header common to all the files being read (i.e. with all read groups merged etc)
   */
  public SAMFileHeader header() {
    return mHeader;
  }

  /**
   * @return the set of SAM files to be read
   */
  public Collection<File> files() {
    return mFiles;
  }

  /**
   * @return the number of threads to use for reading
   */
  public int numThreads() {
    return mNumThreads;
  }

  /**
   * @return the SamFilterParams used for individual record filtering. Fileds involving region restriction
   * (either explicit or via BED files) should not be used, as they have been loaded and resolved in the context already.
   */
  public SamFilterParams filterParams() {
    return mParams;
  }

  /**
   * @return true if the context has range restrictions, false if unrestricted
   */
  public boolean hasRegions() {
    return mReferenceRanges != null;
  }

  /**
   * @return the resolved ranges that the reading should operate over, or null if no ranges are being applied.
   */
  public ReferenceRanges referenceRanges() {
    return mReferenceRanges;
  }

  /**
   * @param sequenceId the SAM sequence id of interest.
   * @return the resolved ranges for the specified template sequence
   * @throws java.lang.NullPointerException if no ranges are being applied
   */
  public RangeList<String> rangeList(int sequenceId) {
    return mReferenceRanges.get(sequenceId);
  }

  /**
   * @param sequenceName the sequence name of interest.
   * @return the resolved ranges for the specified template sequence, or null if no ranges are being applied
   * @throws java.lang.NullPointerException if no ranges are being applied
   */
  public RangeList<String> rangeList(String sequenceName) {
    return mReferenceRanges.get(sequenceName);
  }
}
