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

import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.CRC32;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.io.FileUtils;

/**
 * Helper class for writing sequence data and pointers.
 */
@TestClass(value = {"com.rtg.reader.SequencesWriterTest"})
class NormalSequenceFilePair implements SequenceFilePair {
  private final OutputStream mSeqData;
  private final OutputStream mQualData;
  private final DataOutputStream mPointers;
  private final long mLimit;
  private long mSeqSize = 0;
  private long mQualSize = 0;
  private long mPointerSize = 0;

  private final CRC32 mChecksumSeq;
  private final CRC32 mChecksumQual;
  private final int mPointerEntrySize;


  public NormalSequenceFilePair(OutputStream seqData, OutputStream qualData, DataOutputStream pointers, boolean quality, long limit, CRC32 checksumSeq, CRC32 checksumQual) throws IOException {
    if (limit > Integer.MAX_VALUE) {
      try {
        try {
          seqData.close();
        } finally {
          if (qualData != null) {
            qualData.close();
          }
        }
      } finally {
        pointers.close();
      }
      //TODO handle larger spaces
      throw new IllegalArgumentException("Currently only support int pointers");
    }
    mSeqData =  seqData;
    mQualData = qualData;
    mPointers = pointers;
    mLimit = limit;

    mChecksumSeq = checksumSeq;
    mChecksumQual = quality ? checksumQual : null;
    mPointerEntrySize = 4 + 1 + (quality ? 1 : 0);
  }

  /**
   * Constructor
   * @param dir directory to write files to
   * @param fileNum file number to write
   * @param quality whether to write quality data
   * @param limit largest number of values allowed
   * @param checksumSeq checksum tracker for sequence data
   * @param checksumQual checksum tracker for quality data
   * @throws IOException if an IO error occurs
   */
  public NormalSequenceFilePair(File dir, int fileNum, boolean quality, final long limit, CRC32 checksumSeq, CRC32 checksumQual) throws IOException {
    this(FileUtils.createOutputStream(SdfFileUtils.sequenceDataFile(dir, fileNum), false),
          quality ? FileUtils.createOutputStream(SdfFileUtils.qualityDataFile(dir, fileNum), false) : null,
          new DataOutputStream(FileUtils.createOutputStream(SdfFileUtils.sequencePointerFile(dir, fileNum), false)),
          quality, limit, checksumSeq, checksumQual);

  }

  /**
   * Declares the beginning of a sequence.
   * @return false if operation was not performed due to not being enough space.
   * @throws IOException if an IO error occurs.
   */
  @Override
  public boolean markNextSequence() throws IOException {
    if (mSeqSize < mLimit && mPointerSize + mPointerEntrySize <= mLimit) {
      mPointers.writeByte((int) mChecksumSeq.getValue());
      mChecksumSeq.reset();
      if (mChecksumQual != null) {
        mPointers.writeByte((int) mChecksumQual.getValue());
        mChecksumQual.reset();
      }
      mPointers.writeInt((int) mSeqSize);
      mPointerSize += mPointerEntrySize;
      return true;
    }
    return false;
  }

  @Override
  public void lastSequence() throws IOException {
    mPointers.writeByte((int) mChecksumSeq.getValue());
    mChecksumSeq.reset();
    if (mChecksumQual != null) {
      mPointers.writeByte((int) mChecksumQual.getValue());
      mChecksumQual.reset();
    }
  }


  /**
   * Writes sequence bytes to disk.
   * @param data byte array to write
   * @param offset position in the array to start from (zero based)
   * @param length amount of data to write
   * @return true if succeeded, false if space limit hit.
   * @throws IOException if an I/O error occurs.
   */
  @Override
  public boolean write(byte[] data, int offset, int length) throws IOException {
    if (mSeqSize + length <= mLimit) {
      mSeqData.write(data, offset, length);
      mSeqSize += length;
      mChecksumSeq.update(data, offset, length);
      return true;
    }
    return false;
  }

  /**
   * Writes quality bytes to disk.
   * @param qual byte array to write
   * @param offset position in the array to start from (zero based)
   * @param length amount of data to write
   * @return true if succeeded, false if space limit hit.
   * @throws IOException if an I/O error occurs.
   */
  @Override
  public boolean writeQuality(byte[] qual, int offset, int length) throws IOException {
    if (mQualSize + length <= mLimit) {
      if (length > 0) {
        mQualData.write(qual, offset, length);
        mQualSize += length;
        mChecksumQual.update(qual, offset, length);
      }
      return true;
    }
    return false;
  }

  /**
   * @return original length of data written
   */
  @Override
  public long valuesWritten() {
    return mSeqSize;
  }

  /**
   * @return the number of bytes available to write
   */
  @Override
  public long bytesFree() {
    return mLimit - mSeqSize;
  }

  /**
   * Closes OutputStreams
   * @throws IOException if an I/O error occurs.
   */
  @Override
  public void close() throws IOException {
    mSeqData.close();
    mPointers.close();
    if (mQualData != null) {
      mQualData.close();
    }
  }

}
