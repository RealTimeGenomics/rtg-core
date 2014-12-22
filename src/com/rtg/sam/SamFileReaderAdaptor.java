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

import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.intervals.ReferenceRanges;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

/**
 * Wraps around a <code>SAMFileReader</code> iterator, supplying simple region filtering if needed.
 */
@TestClass("com.rtg.sam.SamFileAndRecordTest")
public class SamFileReaderAdaptor extends AbstractSamRecordIterator {
  private final SAMFileReader mReader;
  private final CloseableIterator<SAMRecord> mIterator;
  private boolean mIsClosed;

  /**
   * Constructor that adapts a regular <code>SAMFileReader</code>, optionally filtering on regions
   * @param reader the reader
   * @param regions regions to filter from output, may be null
   */
  public SamFileReaderAdaptor(SAMFileReader reader, final ReferenceRanges regions) {
    super(reader.getFileHeader());
    mReader = reader;
    SamUtils.logRunId(mReader.getFileHeader());
    if (regions == null || regions.allAvailable()) {
      mIterator = mReader.iterator();
    } else {
      mIterator = new SamRestrictingIterator(mReader.iterator(), regions);
    }
  }

  @Override
  public void close() throws IOException {
    if (!mIsClosed) {
      mIsClosed = true;
      mReader.close();
    }
  }

  @Override
  public boolean hasNext() {
    return mIterator.hasNext();
  }

  @Override
  public SAMRecord next() {
    return mIterator.next();
  }
}
