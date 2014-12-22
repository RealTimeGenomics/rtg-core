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

package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.sam.SamFilter;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Sequence data source for the pre-mapped SAM and BAM file inputs from Illumina
 */
public final class MappedSamBamSequenceDataSource extends SamBamSequenceDataSource {

  private final Map<String, SamSequence> mRecordMap;
  private long mDuplicates = 0;

  private MappedSamBamSequenceDataSource(FileStreamIterator inputs, boolean paired, boolean flattenPaired, SamFilter filter) {
    super(inputs, paired, flattenPaired, filter);
    if (paired) {
      mRecordMap = new HashMap<>();
    } else {
      mRecordMap = null;
    }
  }

  /**
   * Construct a pre-mapped SAM or BAM sequence data source from list of SAM or BAM files
   * @param files list of the SAM or BAM file to use as a sequence data source
   * @param paired true if input will be paired, false otherwise
   * @param flattenPaired if <code>paired</code> is false then this will load both arms into a single SDF
   * @param filter this filter will be applied to the sam records
   * @return SamBamSequenceDataSource the sequence data source for the inputs
   */
  public static MappedSamBamSequenceDataSource fromInputFiles(List<File> files, boolean paired, boolean flattenPaired, SamFilter filter) {
    return new MappedSamBamSequenceDataSource(new FileStreamIterator(files, null), paired, flattenPaired, filter);
  }

  @Override
  protected void checkSortOrder() { }

  @Override
  protected boolean nextRecords() throws IOException {
    if (mPaired) {
      SamSequence rec;
      while ((rec = nextRecord()) != null) {
        checkRecordPaired(rec);
        final SamSequence pair = mRecordMap.remove(rec.getReadName());
        if (pair == null) {
          mRecordMap.put(rec.getReadName(), rec);
        } else {
          placePairedRecord(rec);
          placePairedRecord(pair);
          if (mRecords[0] == null || mRecords[1] == null) {
            final SamSequence r = mRecords[0] == null ? mRecords[1] : mRecords[0];

            if (mDuplicates < 5) {
              Diagnostic.warning("Read " + r.getReadName() + " is duplicated in SAM input.");
              if (mDuplicates == 4) {
                Diagnostic.warning("Subsequent warnings of this type will not be shown.");
              }
            }
            mRecordMap.put(rec.getReadName(), rec);
            mDuplicates++;

            continue;
          }
          return haveNextRecords();
        }
      }
      if (mRecordMap.size() != 0) {
        //throw new NoTalkbackSlimException(mRecordMap.size() + " reads missing a pair when processing paired end SAM input.");
        Diagnostic.warning(mRecordMap.size() + " reads missing a pair when processing paired end SAM input.");
      }
      if (mDuplicates > 0) {
        Diagnostic.warning(mDuplicates + " records ignored as duplicates in input");
      }
      return false;
    } else {
      return super.nextRecords();
    }
  }
}
