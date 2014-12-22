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

import com.rtg.util.io.SeekableStream;

/**
 * Stream manager for reading data, should only be used by DefaultSequencesReader
 */
class SequenceStreamManager extends AbstractStreamManager {

  protected RollingFile mQuality;

  protected boolean mOpenQuality;

  protected final PointerFileHandler mPointerHandler;

  /**
   * Creates a the stream manager
   * @param dir Directory containing sequence data
   * @param numberSequences Number of sequences in directory
   * @throws IOException If an I/O Error occurs
   */
  SequenceStreamManager(final File dir, final long numberSequences, final boolean quality, IndexFile mainIndex, DataFileOpenerFactory openerFactory) throws IOException {
    super(dir, numberSequences, SdfFileUtils.SEQUENCE_INDEX_FILENAME, SdfFileUtils.SEQUENCE_DATA_FILENAME, SdfFileUtils.SEQUENCE_POINTER_FILENAME, mainIndex.dataIndexVersion(), openerFactory.getSequenceOpener());
    mOpenQuality = quality;
    if (mOpenQuality) {
      mQuality = new DataRollingFile(mDir, SdfFileUtils.SEQUENCE_QUALITY_DATA_FILENAME, mIndex.numberEntries(), mIndex, openerFactory.getQualityOpener());
    }
    mPointerHandler = PointerFileHandler.getHandler(mainIndex, PointerFileHandler.SEQUENCE_POINTER);
  }

  @Override
  protected void seekImpl(final long seqNum) throws IOException {
    final long pointerpos = seqNum - mCurrentLower;
    mPointerHandler.initialisePosition(mIndexedSequenceFileNumber, pointerpos, mPointers, mIndex);
    if (mPointerHandler.seqLength() > Integer.MAX_VALUE || mPointerHandler.seqLength() < 0) {
      //we only allow single sequences up to 2gb
      throw new CorruptSdfException();
    }
    mDataLength = mPointerHandler.seqLength();
    mData.randomAccessFile().seek(mPointerHandler.seqPosition());
    if (mOpenQuality) {
      mQuality.randomAccessFile().seek(mPointerHandler.seqPosition());
    }
  }

  /**
   */
  @Override
  protected void openFiles() throws IOException {
    super.openFiles();
    if (mOpenQuality) {
      openQualityFile(mIndexedSequenceFileNumber);
    }
  }

  protected void openQualityFile(final int fileno) throws IOException {
    if (!mQuality.openDataFile(fileno)) {
      throw new CorruptSdfException("Expected file missing");
    }
  }

  @Override
  protected void ensureFiles() throws IOException {
    super.ensureFiles();
    if (mOpenQuality && mQuality.currentFileNo() != mIndexedSequenceFileNumber) {
      openQualityFile(mIndexedSequenceFileNumber);
    }
  }

  @Override
  public void close() throws IOException {
    super.close();
    if (mOpenQuality && mQuality != null) {
      mQuality.close();
    }
  }

  int read(byte[] dataOut, int start, int length, RollingFile src) throws IOException {
    if (length == 0) {
      return 0;
    }
    final int fullLength = (int) getDataLength();
    if ((start + length) > fullLength) {
      throw new IllegalArgumentException("Requested data not a subset of sequence data.");
    }
    if (length > dataOut.length) {
      throw new IllegalArgumentException("Array too small got: " + dataOut.length + " required: " + length);
    }
    SeekableStream raf = src.randomAccessFile();
    int count = 0;
    while ((start - count) >= (int) (raf.length() - raf.getPosition())) {
      count += (int) (raf.length() - raf.getPosition());
      if (!src.rollFile()) {
        throw new CorruptSdfException("expected SDF data missing");
      }
      raf = src.randomAccessFile();
    }
    raf.seek(raf.getPosition() + (start - count));
    count = 0;
    int last;
    while (count < length) {
      last = raf.read(dataOut, count, length - count);
      final int eof = -1;
      if (last == eof) {
        //src.rollFile();
        if (!src.rollFile()) {
          throw new CorruptSdfException("expected SDF data missing");
        }
        raf = src.randomAccessFile();
        last = raf.read(dataOut, count, length - count);
      }
      if (last >= 0) {
        count += last;
      }
    }
    return length;
  }

  int readData(byte[] dataOut, int start, int length) throws IOException {
    return read(dataOut, start, length, mData);
  }

  int readQuality(byte[] dataOut, int start, int length) throws IOException {
    return read(dataOut, start, length, mQuality);
  }

  byte sequenceChecksum() {
    return mPointerHandler.checksum();
  }
  byte qualityChecksum() {
    return mPointerHandler.qualityChecksum();
  }
}
