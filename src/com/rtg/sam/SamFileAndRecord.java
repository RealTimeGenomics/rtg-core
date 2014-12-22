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

import com.rtg.launcher.GlobalFlags;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * Adds some support to SkipInvalidRecordsIterator to allow externally coordinated reading of multiple SAM files.
 */

public class SamFileAndRecord extends SkipInvalidRecordsIterator implements Comparable<SamFileAndRecord>, Integrity {

  private final int mId;
  private final String mPath;

  // work around to allow user to force file loading when sam header is malformed
  private static final boolean IGNORE_HEADER_SORTORDER = GlobalFlags.isSet(GlobalFlags.SAM_IGNORE_SORT_ORDER_FLAG);

  protected SamFileAndRecord(String path, int id, RecordIterator<SAMRecord> adaptor) throws IOException {
    super(path, adaptor, false);
    mId = id;
    mPath = path;
    validateParams(path);
  }

  /**
   * @param samFile file to be read.
   * @param id unique identifier within one instance of a <code>SAMMultifileIterator</code>.
   * @throws IOException if an IO Error occurs
   */
  public SamFileAndRecord(final File samFile, final int id) throws IOException {
    super(samFile);
    mPath = samFile.getPath();
    mId = id;
    validateParams(samFile.getPath());
  }

  private void validateParams(String path) throws IOException {
    if (mId < 0) {
      close();
      throw new IllegalArgumentException("id less than zero supplied: " + mId);
    }
    if (header().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
      if (!IGNORE_HEADER_SORTORDER) {
        throw new IOException("\"" + path + "\" is not sorted in coordinate order.");
      }
    }
  }

  /**
   * Get the name of the underlying file.
   * @return the name of the underlying file.
   */
  public String name() {
    return mPath;
  }

  @Override
  public boolean equals(Object obj) {
    // to keep findbugs quiet
    return super.equals(obj);
  }

  @Override
  public int hashCode() {
    // to keep findbugs quiet
    return super.hashCode();
  }

  @Override
  public int compareTo(final SamFileAndRecord that) {
    final int samCompare = SamCompareUtils.compareSamRecords(this.mRecord, that.mRecord);
    if (samCompare != 0) {
      return samCompare;
    }
    final int thisId = this.mId;
    final int thatId = that.mId;
    //System.err.println("id  this=" + thisId + " that=" + thatId);
    if (thisId < thatId) {
      return -1;
    }
    if (thisId > thatId) {
      return +1;
    }
    assert this == that;
    return 0;
  }

  String orderToString() {
    return ("index=" + mRecord.getReferenceIndex()) + " start=" + mRecord.getAlignmentStart() + " id=" + mId;
  }

  @Override
  public String toString() {
    return "SamFileAndRecord:" + mId + " line=" + mRecordCount + " file=" + mPath + " record=" + mRecord.getSAMString().trim();
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mId >= 0);

    return true;
  }

  @Override
  public boolean globalIntegrity() {
    return integrity();
  }

}
