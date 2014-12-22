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
import java.io.RandomAccessFile;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.ByteArrayIOUtils;

/**
 * Handles reading of a data file index
 */
@TestClass("com.rtg.reader.SequenceStreamManagerTest")
public abstract class DataFileIndex {

  static final long SEQ_COUNT_ONLY_VERSION = 1L;
  static final long DATASIZE_VERSION = 2L;

  DataFileIndex() {
  }

  abstract long getTotalNumberSequences();
  abstract long numberSequences(int fileIndex);
  abstract long dataSize(int fileIndex);
  abstract int numberEntries();


  static DataFileIndex loadSequenceDataFileIndex(long version, File dir) throws IOException {
    return loadDataFileIndex(version, SdfFileUtils.sequenceIndexFile(dir), SdfFileUtils.SEQUENCE_DATA_FILENAME);
  }

  static DataFileIndex loadLabelDataFileIndex(long version, File dir) throws IOException {
    return loadDataFileIndex(version, SdfFileUtils.labelIndexFile(dir), SdfFileUtils.LABEL_DATA_FILENAME);
  }

  static DataFileIndex loadLabelSuffixDataFileIndex(long version, File dir) throws IOException {
    return loadDataFileIndex(version, SdfFileUtils.labelSuffixIndexFile(dir), SdfFileUtils.LABEL_SUFFIX_DATA_FILENAME);
  }

  static DataFileIndex loadDataFileIndex(long version, File dataIndexFile, String dataFilePrefix) throws IOException {
    if (version == SEQ_COUNT_ONLY_VERSION) {
      return new DataFileIndexVersion1(dataIndexFile, dataFilePrefix);
    } else if (version == DATASIZE_VERSION) {
      return new DataFileIndexVersion2(dataIndexFile);
    } else {
      throw new NoTalkbackSlimException("Unsupported SDF version: " + dataIndexFile.getParent()); //should be unpossible
    }
  }

  @TestClass("com.rtg.reader.SequenceStreamManagerTest")
  private static class DataFileIndexVersion2 extends DataFileIndex {
    private final long[] mIndex;
    private final long mTotal;
    private final int mNumberEntries;

    public DataFileIndexVersion2(File dataIndexFile) throws IOException {
      if (dataIndexFile.exists()) {
        try (RandomAccessFile indexRAF = new RandomAccessFile(dataIndexFile, "r")) {
          //TODO handle longer indexs
          //current maximum by this implementation would be:
          //Integer.MAX_VALUE / bytes per entry = max_entries
          //2147483647        / 16              = 134217727
          //max_entries * data_size_per_entry   = max_sdf_size
          //134217727   * 2147483648 (default)  = 288230374004228096 (262143 TB) of uncompressed 1 byte per nucleotide data
          final byte[] tempData = new byte[(int) indexRAF.length()];
          indexRAF.readFully(tempData);
          try {
            mIndex = ByteArrayIOUtils.convertToLongArray(tempData); //new long[(int) (mIndexRAF.length() / 8)];
          } catch (final IllegalArgumentException e) {
            throw new CorruptSdfException(dataIndexFile.getParentFile());
          }
          long count = 0;
          for (int i = 0; i < mIndex.length; i += 2) {
            //mIndex[i] = mIndexRAF.readLong();
            if (mIndex[i] < 0) {
              throw new CorruptSdfException(dataIndexFile.getParentFile());
            }
            count += mIndex[i];
          }
          mTotal = count;
          mNumberEntries = mIndex.length / 2;
        }
      } else {
        mIndex = new long[0];
        mTotal = 0;
        mNumberEntries = 0;
      }
    }
    @Override
    long dataSize(int fileIndex) {
      return mIndex[fileIndex * 2 + 1];
    }
    @Override
    long getTotalNumberSequences() {
      return mTotal;
    }
    @Override
    int numberEntries() {
      return mNumberEntries;
    }
    @Override
    long numberSequences(int fileIndex) {
      return mIndex[fileIndex * 2];
    }
  }

  @TestClass("com.rtg.reader.SequenceStreamManagerTest")
  private static class DataFileIndexVersion1 extends DataFileIndex {

    private final File mDir;
    private final long[] mIndex;
    private final long mTotal;
    private final String mDataFilePrefix;

    public DataFileIndexVersion1(File dataIndexFile, String dataFilePrefix) throws IOException {
      mDataFilePrefix = dataFilePrefix;
      mDir = dataIndexFile.getParentFile();
      if (dataIndexFile.exists()) {
        try (RandomAccessFile indexRAF = new RandomAccessFile(dataIndexFile, "r")) {
          final byte[] tempData = new byte[(int) indexRAF.length()];
          indexRAF.readFully(tempData);
          try {
            mIndex = ByteArrayIOUtils.convertToLongArray(tempData); //new long[(int) (mIndexRAF.length() / 8)];
          } catch (final IllegalArgumentException e) {
            throw new CorruptSdfException(mDir);
          }
          long count = 0;
          for (final long i : mIndex) {
            if (i < 0) {
              throw new CorruptSdfException(mDir);
            }
            count += i;
          }
          mTotal = count;
        }
      } else {
        mIndex = new long[0];
        mTotal = 0;
      }

    }
    @Override
    long dataSize(int fileIndex) {
      return new File(mDir, mDataFilePrefix + fileIndex).length();
    }
    @Override
    long numberSequences(int fileIndex) {
      return mIndex[fileIndex];
    }
    @Override
    long getTotalNumberSequences() {
      return mTotal;
    }
    @Override
    int numberEntries() {
      return mIndex.length;
    }
  }
}
