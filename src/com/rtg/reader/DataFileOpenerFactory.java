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
import java.io.InputStream;

import com.rtg.mode.SequenceType;
import com.rtg.util.bytecompression.CompressedByteArray;
import com.rtg.util.io.BufferedRandomAccessFile;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.RandomAccessFileStream;
import com.rtg.util.io.SeekableStream;

/**
 * Class for opening variety of SDF data files
 */
public class DataFileOpenerFactory {

  private final DataFileOpener mLabel;
  private final DataFileOpener mSequence;
  private final DataFileOpener mQuality;

  /**
   * Create a factory for providing {@link DataFileOpener}'s for given parameters
   * @param sequenceEncoding the sequence encoding (obtained from {@link com.rtg.reader.IndexFile})
   * @param qualityEncoding the quality encoding (obtained from {@link com.rtg.reader.IndexFile})
   * @param type the sequence type
   */
  public DataFileOpenerFactory(byte sequenceEncoding, byte qualityEncoding, SequenceType type) {
    mLabel = new NormalOpener();
    mSequence = IndexFile.SEQUENCE_ENCODING_COMPRESSED == sequenceEncoding ? new BitwiseOpener(CompressedByteArray.minBits(type.numberCodes())) : new NormalOpener();
    mQuality = IndexFile.QUALITY_ENCODING_COMPRESSED == qualityEncoding ? new CompressedOpener(CompressedMemorySequencesReader.MAX_QUAL_VALUE) : new NormalOpener();
  }

  /**
   * @return opener for label data
   */
  public DataFileOpener getLabelOpener() {
    return mLabel;
  }

  /**
   * @return opener for sequence data
   */
  public DataFileOpener getSequenceOpener() {
    return mSequence;
  }

  /**
   * @return opener for quality data
   */
  public DataFileOpener getQualityOpener() {
    return mQuality;
  }

  private static class NormalOpener implements DataFileOpener {
    @Override
    public InputStream open(File filename, long length) throws IOException {
      return FileUtils.createFileInputStream(filename, false);
    }
    @Override
    public SeekableStream openRandomAccess(File filename, long length) throws IOException {
      return new RandomAccessFileStream(new BufferedRandomAccessFile(filename, "r"));
    }
  }
  private static class BitwiseOpener implements DataFileOpener {
    private final int mBits;
    public BitwiseOpener(int bits) {
      mBits = bits;
    }
    @Override
    public InputStream open(File filename, long length) throws IOException {
      return new FileBitwiseInputStream(filename, mBits, length, false);
    }
    @Override
    public SeekableStream openRandomAccess(File filename, long length) throws IOException {
      return new FileBitwiseInputStream(filename, mBits, length, true);
    }
  }

  private static class CompressedOpener implements DataFileOpener {
    private final int mRange;
    public CompressedOpener(int range) {
      mRange = range;
    }
    @Override
    public InputStream open(File filename, long length) throws IOException {
      return new FileCompressedInputStream(filename, mRange, length, false);
    }
    @Override
    public SeekableStream openRandomAccess(File filename, long length) throws IOException {
      return new FileCompressedInputStream(filename, mRange, length, true);
    }
  }
}
