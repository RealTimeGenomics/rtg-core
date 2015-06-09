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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.tabix.LocusIndex;
import com.rtg.tabix.TabixIndexReader;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.VirtualOffsets;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.ClosedFileInputStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Utility methods for getting a <code>SAMRecord</code> iterator that is backed by a {@link ClosedFileInputStream}. Note there is no
 * way to reset the reader. If you want to create a new one you must instantiate a new instance of this class.
 */
public final class SamClosedFileReader extends AbstractSamRecordIterator {

  private static final boolean MULTI_CHUNKS = true;

  private final File mFile;
  private final boolean mIsBam;
  private final ReferenceRanges<String> mRegions;
  private final CloseableIterator<SAMRecord> mIterator;
  private final ClosedFileInputStream mStream;
  private boolean mIsClosed;

  /**
   * @param file SAM or BAM file, if region is specified an index must also be present with the appropriate relative file path
   * @param regions regions the file should select, null for whole file
   * @param header header that should be used for SAM records may not be null
   * @throws IOException if an IO error occurs
   */
  public SamClosedFileReader(File file, ReferenceRanges<String> regions, SAMFileHeader header) throws IOException {
    super(header);
    SamUtils.logRunId(header);
    mIsBam = SamUtils.isBAMFile(file);
    mFile = file;
    mStream = new ClosedFileInputStream(mFile);
    mRegions = regions;
    mIterator = obtainIteratorInternal();
  }

  @Override
  public boolean hasNext() {
    return mIterator.hasNext();
  }

  @Override
  public SAMRecord next() {
    return mIterator.next();
  }

  @Override
  public void close() throws IOException {
    if (!mIsClosed) {
      mIsClosed = true;
      mIterator.close();
    }
  }


  private CloseableIterator<SAMRecord> primaryIterator(BlockCompressedInputStream bcis) throws IOException {
    return SamUtils.makeSamReader(bcis, mHeader, mIsBam ? SamReader.Type.BAM_TYPE : SamReader.Type.SAM_TYPE).iterator();
  }

  /**
   * Creates a closed file input stream backed {@link SAMRecord} iterator over the given region.
   * @return the iterator
   * @throws IOException if an IO error occurs
   */
  private CloseableIterator<SAMRecord> obtainIteratorInternal() throws IOException {

    // Handle easy case of no restriction
    if (mRegions == null || mRegions.allAvailable()) {
      if (!mIsBam) {
        return SamUtils.makeSamReader(mStream, mHeader).iterator(); // htsjdk will decide whether decompression is required
      } else {
        return primaryIterator(new BlockCompressedInputStream(mStream));
      }
    }

    final LocusIndex index;
    if (mIsBam) {
      File indexFileName = BamIndexer.indexFileName(mFile);
      if (!indexFileName.exists()) {
        indexFileName = BamIndexer.secondaryIndexFileName(mFile);
        if (!indexFileName.exists()) {
          throw new NoTalkbackSlimException("File " + mFile.getPath() + " is not indexed");
        }
      }
      index = new BamIndexReader(indexFileName, mHeader.getSequenceDictionary());
    } else {
      final File indexFileName = TabixIndexer.indexFileName(mFile);
      if (!TabixIndexer.isBlockCompressed(mFile) || !indexFileName.exists()) {
        throw new NoTalkbackSlimException("File " + mFile.getPath() + " is not indexed");
      }
      index = new TabixIndexReader(indexFileName);
    }

    final VirtualOffsets filePointers = index.getFilePointers(mRegions);

    if (filePointers == null) {
      return SamUtils.makeSamReader(new ByteArrayInputStream(new byte[0]), mHeader, SamReader.Type.SAM_TYPE).iterator();
    }

    //Diagnostic.developerLog("Using virtual offsets for file: " + mFile.toString() + "\t" + filePointers);

    /*
     * NOTE: now we make sure we handle skipping records that are within the indexed region but not within given restriction,
     * by wrapping in a restricting iterator.
     */
    final BlockCompressedInputStream cfis = new BlockCompressedInputStream(mStream);
    if (MULTI_CHUNKS) {
      return new SamMultiRestrictingIterator(cfis, filePointers, mHeader, mIsBam, mFile.toString());
    } else {
      // This should be correct but will read all data from start of first region through to end of last region
      cfis.seek(filePointers.start(0));
      return new SamRestrictingIterator(primaryIterator(cfis), mRegions);
    }
  }
}
