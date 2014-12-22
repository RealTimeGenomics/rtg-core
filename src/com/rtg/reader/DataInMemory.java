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
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

import com.rtg.mode.SequenceType;
import com.rtg.util.bytecompression.BitwiseByteArray;
import com.rtg.util.bytecompression.CompressedByteArray;

/**
 * Helper for <code>CompressedMemorySequencesReader2</code> does most of the heavy lifting.
 * Loads compressed versions of SDF directly from disk into memory without having to decompress and recompress.
 *
 */
public final class DataInMemory {
  private final int[][] mPointers;
  private final byte[][] mSequenceChecksums;
  private final byte[][] mQualityChecksums;
  private final BitwiseByteArray[] mSequenceData;
  private CompressedByteArray[] mQualityData;
  private final PointerFileLookup mFileNoLookup;

  private final QualityLoader mQualLoader;

  //current sequence
  private int mSeqId = -1;
  private int mStartFileNo;
  private int mStartIndex;
  private int mSeqLength;

  private DataInMemory(int[][] pointers, byte[][] sequenceChecksums, byte[][] qualityChecksums, BitwiseByteArray[] sequenceData, QualityLoader loader) {
    mPointers = pointers;
    mSequenceChecksums = sequenceChecksums;
    mQualityChecksums = qualityChecksums;
    mSequenceData = sequenceData;
    mFileNoLookup = PointerFileLookup.generateLookup(mPointers);
    mQualLoader = loader;
  }

  void initQuality() throws IOException {
    if (mQualityData == null) {
      //mQualLoader is shared between all threads, this DataInMemory instance (and its reference to mQualityData) should only be accessed by one thread ever
      //therefore getting/setting the mQualityData reference doesn't need to be synchronized. The synchronization for initialization is
      //performed within the quality loader. Resulting in this instance having a reference to the shared mQualityData that should only ever be used for read operations.
      mQualityData = mQualLoader.loadQualityData();
    }
  }

  private void seq(int seqId) {
    if (mSeqId != seqId) {
      mSeqId = seqId;
      mStartFileNo = mFileNoLookup.lookup(seqId);
      mStartIndex = seqId - mFileNoLookup.startSeq(mStartFileNo);
      mSeqLength = -1;
    }
  }

  /**
   * @param seqId sequence id
   * @return the length of given sequence
   */
  public int length(int seqId) {
    seq(seqId);
    if (mSeqLength == -1) {
      int tot = mPointers[mStartFileNo][mStartIndex + 1] - mPointers[mStartFileNo][mStartIndex];
      if (mStartIndex + 1 == mPointers[mStartFileNo].length - 1) {
        int fileNo = mStartFileNo + 1;
        while (fileNo < mPointers.length && mPointers[fileNo].length == 1) {
          tot += mPointers[fileNo][0];
          fileNo++;
        }
        if (fileNo < mPointers.length) {
          tot += mPointers[fileNo][0];
        }
      }
      mSeqLength = tot;
    }
    return mSeqLength;
  }

  /**
   * Read sequence data
   * @param seq the sequence id
   * @param dest the destination array for data
   * @param start the start position within the sequence to read
   * @param length the maximum amount of data to read from the sequence.
   * @return the amount of values read
   */
  public int readSequence(int seq, byte[] dest, int start, int length) {
    seq(seq);
    if (mPointers[mStartFileNo][mStartIndex] == mPointers[mStartFileNo][mStartIndex + 1]) {
      // Return 0 if the sequence has 0 length
      return 0;
    }
    int arrOffset = mPointers[mStartFileNo][mStartIndex] + start;
    int fileNo = mStartFileNo;
    int index = mStartIndex + 1;
    while (arrOffset >= mPointers[fileNo][index]) {
      arrOffset -= mPointers[fileNo][index];
      fileNo++;
      index = 0;
    }
    int lengthToGet = Math.min(mPointers[fileNo][index] - arrOffset, length);
    int destOffset = 0;
    mSequenceData[fileNo].get(dest, arrOffset, 0, lengthToGet);
    destOffset += lengthToGet;
    if (destOffset < length && index == mPointers[fileNo].length - 1) {
      fileNo++;
      while (destOffset < length && fileNo < mPointers.length && mPointers[fileNo].length == 1) {
        lengthToGet = Math.min(mPointers[fileNo][0], length - destOffset);
        mSequenceData[fileNo].get(dest, 0, destOffset, lengthToGet);
        destOffset += lengthToGet;
        fileNo++;
      }
      if (destOffset < length && fileNo < mPointers.length) {
        lengthToGet = Math.min(mPointers[fileNo][0], length - destOffset);
        mSequenceData[fileNo].get(dest, 0, destOffset, lengthToGet);
        destOffset += lengthToGet;
      }
    }
    return destOffset;
  }

  /**
   * Read quality data
   * @param seq the sequence id
   * @param dest the destination array for data
   * @param start the start position within the sequence to read
   * @param length the maximum amount of quality data to read from the sequence.
   * @return the amount of values read
   * @throws IOException when an IO error occurs
   */
  public int readQuality(int seq, byte[] dest, int start, int length) throws IOException {
    seq(seq);
    if (mPointers[mStartFileNo][mStartIndex] == mPointers[mStartFileNo][mStartIndex + 1]) {
      // Return 0 if the sequence has 0 length
      return 0;
    }
    initQuality();
    int arrOffset = mPointers[mStartFileNo][mStartIndex] + start;
    int fileNo = mStartFileNo;
    int index = mStartIndex + 1;
    while (arrOffset >= mPointers[fileNo][index]) {
      arrOffset -= mPointers[fileNo][index];
      fileNo++;
      index = 0;
    }
    int lengthToGet = Math.min(mPointers[fileNo][index] - arrOffset, length);
    int destOffset = 0;
    mQualityData[fileNo].get(dest, arrOffset, 0, lengthToGet);
    destOffset += lengthToGet;
    if (destOffset < length && index == mPointers[fileNo].length - 1) {
      fileNo++;
      while (destOffset < length && fileNo < mPointers.length && mPointers[fileNo].length == 1) {
        lengthToGet = Math.min(mPointers[fileNo][0], length - destOffset);
        mQualityData[fileNo].get(dest, 0, destOffset, lengthToGet);
        destOffset += lengthToGet;
        fileNo++;
      }
      if (destOffset < length && fileNo < mPointers.length) {
        lengthToGet = Math.min(mPointers[fileNo][0], length - destOffset);
        mQualityData[fileNo].get(dest, 0, destOffset, lengthToGet);
        destOffset += lengthToGet;
      }
    }
    return destOffset;
  }

  /**
   * return total length between given sequences
   * @param start sequence id (inclusive)
   * @param end sequence id (exclusive)
   * @return total length
   */
  public long lengthBetween(int start, int end) {
    final int fileStart = mFileNoLookup.lookup(start);
    final int fileEnd = mFileNoLookup.lookup(end);
    final int startIndex = start - mFileNoLookup.startSeq(fileStart);
    final int endIndex = end - mFileNoLookup.startSeq(fileEnd);
    long retLength = -mPointers[fileStart][startIndex];
    for (int i = fileStart; i < fileEnd; i++) {
      retLength += mPointers[i][mPointers[i].length - 1];
    }
    retLength += mPointers[fileEnd][endIndex];
    return retLength;
  }

  /**
   * return array of sequence lengths
   * @param start sequence id (inclusive)
   * @param end sequence id (exclusive)
   * @return array of lengths
   */
  public int[] sequenceLengths(int start, int end) {
    final int fileStart = mFileNoLookup.lookup(start);
    final int startIndex = start - mFileNoLookup.startSeq(fileStart);
    final int[] ret = new int[end - start];
    int fileNo = fileStart;
    int pointerIndex = startIndex;
    for (int retIndex = 0; retIndex < end - start; retIndex++) {
      ret[retIndex] = -mPointers[fileNo][pointerIndex];
      pointerIndex++;
      while (fileNo < mPointers.length && pointerIndex == mPointers[fileNo].length - 1) {
        ret[retIndex] += mPointers[fileNo][pointerIndex];
        fileNo++;
        pointerIndex = 0;
      }
      if (fileNo < mPointers.length && pointerIndex < mPointers[fileNo].length) {
        ret[retIndex] += mPointers[fileNo][pointerIndex];
      }
    }
    return ret;
  }

  private static int numberEntries(DataFileIndex seqIndex, long start, long end) {
    int numberEntries = 0;
    long runningTotal = 0;
    for (int i = 0; i < seqIndex.numberEntries(); i++) {
      final long addSeqs = seqIndex.numberSequences(i);
      if (start < runningTotal + addSeqs && end >= runningTotal) {
        numberEntries++;
      }
      runningTotal += addSeqs;
    }
    return numberEntries;
  }

  /**
   * create a shallow copy
   * @return the copy
   */
  public DataInMemory copy() {
    return new DataInMemory(mPointers, mSequenceChecksums, mQualityChecksums, mSequenceData, mQualLoader);
  }

    /**
   * create data in memory
   * @param directory directory of SDF
   * @param index index file
   * @param seqIndex sequence index
   * @param start sequence id (inclusive)
   * @param end sequence id (exclusive)
   * @return the data helper
   * @throws IOException IO error occurs
   */
  public static DataInMemory loadDelayQuality(File directory, IndexFile index, DataFileIndex seqIndex, long start, long end) throws IOException {
    final PointerLoader point = new PointerLoader(directory, index, seqIndex, start, end);
    final SequenceLoader seq;
    final QualityLoader qual;
    if (start < end) {
      point.loadPointers();
      seq = new SequenceLoader(point);
      qual = new QualityLoader(point);
      seq.loadSequenceData();
    } else {
      seq = new SequenceLoader(point);
      qual = new QualityLoader(point);
    }
    return new DataInMemory(point.mPointers, point.mSequenceChecksums, point.mQualityChecksums, seq.mSequenceData, qual);
  }

  static final class PointerLoader {
    private final File mDir;
    private final IndexFile mIndexFile;
    private final boolean mHasQuality;
    private final boolean mHasChecksums;
    private final DataFileIndex mSeqIndex;
    private final long mStart;
    private final long mEnd;
    final int[][] mPointers;
    private final byte[][] mSequenceChecksums;
    private final byte[][] mQualityChecksums;
    private final int mPointerEntrySize;
    private final int mNumberEntries;
    private int mFirstFilePointerAdj = 0;

    PointerLoader(File directory, IndexFile index, DataFileIndex seqIndex, long start, long end) {
      mDir = directory;
      mHasQuality = index.hasQuality();
      mIndexFile = index;
      mSeqIndex = seqIndex;
      mStart = start;
      mEnd = end;
      mNumberEntries = numberEntries(seqIndex, start, end);
      mPointers = new int[mNumberEntries][];
      final int pointerSize = 4;
      final int checksumSize;
      if (index.getVersion() < IndexFile.PER_SEQUENCE_CHECKSUM_VERSION) {
        checksumSize = 0;
        mHasChecksums = false;
      } else {
        mHasChecksums = true;
        checksumSize = mHasQuality ? 2 : 1;
      }
      mSequenceChecksums = mHasChecksums ? new byte[mNumberEntries][] : null;
      mQualityChecksums = mHasChecksums && mHasQuality ? new byte[mNumberEntries][] : null;
      mPointerEntrySize = checksumSize + pointerSize;
    }

    public void loadPointers() throws IOException {
      final ByteBuffer buf = ByteBuffer.allocate(64 * 1024).order(ByteOrder.BIG_ENDIAN);
      long fileStart = 0;
      long fileEnd = 0;
      int i = 0;
      for (int fn = 0; fn < mSeqIndex.numberEntries(); fn++) {
        fileEnd += mSeqIndex.numberSequences(fn);
        if (mStart < fileEnd && mEnd >= fileStart) {
          int j = 0;
          final File f = SdfFileUtils.sequencePointerFile(mDir, fn);
          mPointers[i] = new int[(int) (Math.min(fileEnd, mEnd) - Math.max(fileStart, mStart)) + 1]; //should never have more than 2b entries in a file
          if (mHasChecksums) {
            mSequenceChecksums[i] = new byte[mPointers[i].length];
            if (mHasQuality) {
              mQualityChecksums[i] = new byte[mPointers[i].length];
            }
          }
          try (FileInputStream fis = new FileInputStream(f)) {
            buf.clear();
            final FileChannel channel = fis.getChannel();
            try {
              if (mStart > fileStart) {
                channel.position((mStart - fileStart) * mPointerEntrySize);
              }
              chanread:
              while (channel.read(buf) != -1) {
                buf.flip();
                while (buf.remaining() >= mPointerEntrySize) {
                  if (mHasChecksums) {
                    mSequenceChecksums[i][j] = buf.get();
                    if (mHasQuality) {
                      mQualityChecksums[i][j] = buf.get();
                    }
                  }
                  mPointers[i][j] = buf.getInt();
                  j++;
                  if (j >= mPointers[i].length) {
                    break chanread;
                  }
                }
                buf.compact();
              }
              assert j >= mPointers[i].length - 1;
              if (j == mPointers[i].length - 1) {
                mPointers[i][j] = (int) mSeqIndex.dataSize(fn); //I hates teh cast. But this should be fine due to our internal rules about file sizes.
              }
            } finally {
              channel.close();
            }
          }
          i++;
        }
        fileStart = fileEnd;
      }
      readjustPointers();
    }

    public void readjustPointers() {
      if (mPointers.length > 0) {
        if (mPointers[0][0] > 0) {
          final int adj = mPointers[0][0];
          mFirstFilePointerAdj = adj;
          for (int j = 0; j < mPointers[0].length; j++) {
            mPointers[0][j] -= adj;
          }
        }
      }
    }
  }
  abstract static class DataLoader {
    final int[][] mPointers;
    final int mNumberEntries;
    final IndexFile mIndexFile;
    final DataFileIndex mDataFileIndex;
    final int mFirstFilePointerAdjust;
    final long mStart;
    final long mEnd;
    final File mDir;

    DataLoader(PointerLoader pointers) {
      mPointers = pointers.mPointers;
      mNumberEntries = pointers.mNumberEntries;
      mIndexFile = pointers.mIndexFile;
      mDataFileIndex = pointers.mSeqIndex;
      mFirstFilePointerAdjust = pointers.mFirstFilePointerAdj;
      mStart = pointers.mStart;
      mEnd = pointers.mEnd;
      mDir = pointers.mDir;
    }
  }


  static final class QualityLoader extends  DataLoader {
    private CompressedByteArray[] mQualityData;

    QualityLoader(PointerLoader loader) {
      super(loader);
    }

    public synchronized CompressedByteArray[] loadQualityData() throws IOException {
      if (mIndexFile.hasQuality() && mQualityData == null) {
        final CompressedByteArray[] ret = new CompressedByteArray[mNumberEntries];
        long fileStart = 0;
        long fileEnd = 0;
        int i = 0;
        for (int fn = 0; fn < mDataFileIndex.numberEntries(); fn++) {
          fileEnd += mDataFileIndex.numberSequences(fn);
          if (mStart < fileEnd && mEnd >= fileStart) {
            final int startVal;
            final int endVal;
            if (i == 0) {
              startVal = mPointers[i][0] + mFirstFilePointerAdjust;
              endVal = mPointers[i][mPointers[i].length - 1] + mFirstFilePointerAdjust;
            } else {
              startVal = 0;
              endVal = mPointers[i][mPointers[i].length - 1] + mFirstFilePointerAdjust;
            }
            ret[i] = CompressedByteArray.loadCompressed(SdfFileUtils.qualityDataFile(mDir, fn), startVal, endVal, CompressedMemorySequencesReader.MAX_QUAL_VALUE);
            i++;
          }
          fileStart = fileEnd;
        }
        mQualityData = ret;
      }
      return mQualityData;
    }
  }

  static final class SequenceLoader extends DataLoader {
    final BitwiseByteArray[] mSequenceData;
    SequenceLoader(PointerLoader loader) {
      super(loader);
      mSequenceData = new BitwiseByteArray[mNumberEntries];
    }

    public void loadSequenceData() throws IOException {
      final SequenceType type = SequenceType.values()[mIndexFile.getSequenceType()];
      final int range = type.numberKnownCodes() + type.firstValid();
      long fileStart = 0;
      long fileEnd = 0;
      int i = 0;
      for (int fn = 0; fn < mDataFileIndex.numberEntries(); fn++) {
        fileEnd += mDataFileIndex.numberSequences(fn);
        if (mStart < fileEnd && mEnd >= fileStart) {
          final int startVal;
          final int endVal;
          if (i == 0) {
              startVal = mPointers[i][0] + mFirstFilePointerAdjust;
              endVal = mPointers[i][mPointers[i].length - 1] + mFirstFilePointerAdjust;
            } else {
              startVal = 0;
              endVal = mPointers[i][mPointers[i].length - 1] + mFirstFilePointerAdjust;
            }
          mSequenceData[i] = BitwiseByteArray.loadBitwise(SdfFileUtils.sequenceDataFile(mDir, fn), startVal, endVal, range);
          i++;
        }
        fileStart = fileEnd;
      }
    }
  }
}
