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


import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.regex.Pattern;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.SequenceType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.gzip.GzipUtils;
import com.rtg.util.io.AsynchInputStream;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.bzip2.CBZip2InputStream;

/**
 * Read read data from Complete Genomics ASCII 2.0 format.
 */
public class TsvSequenceDataSource implements CgSequenceDataSource {

  private static final int READ_FIELD = 1;
  private static final int QUALITY_FIELD = 2;

  private InputStream mStream;
  private BufferedReader mReader;
  private final Integer mMaximumNs;
  private final File mReads;
  private int mMaximumNsCount = 0;
  private int mSkippedNReads = 0;
  private long mSkippedNResidues = 0;
  private long mMaxLength = 0;
  private long mMinLength = Long.MAX_VALUE;

  private long mSeqNo = 0;
  private byte[] mLeftRead;
  private byte[] mRightRead;
  private byte[] mLeftQuality;
  private byte[] mRightQuality;

  private int mCurrentLength = 0;

  private static final Pattern PATTERN_TAB_SPLIT = Pattern.compile("\t");

  private final byte[] mFastaSymbolLookupTable = new DNAFastaSymbolTable().getAsciiToOrdinalTable();

  /**
   * Constructs the data source
   * @param reads file containing read data
   * @param maximumNs maximum number of Ns allowed in either side for a read. Null for no filtering.
   */
  public TsvSequenceDataSource(final File reads, final Integer maximumNs) {
    mMaximumNs = maximumNs;
    mReads = reads;
  }

  @Override
  @SuppressWarnings("try")
  public void close() throws IOException {
    try (BufferedReader ignored = mReader) {
      mReader = null;
    } finally {
      try (InputStream ignored = mStream) {
        mStream = null;
      }
    }
  }

  @Override
  public SequenceType type() {
    return SequenceType.DNA;
  }

  @Override
  public boolean hasQualityData() {
    return true;
  }

  private boolean acceptableNs(final String s) {
    if (mMaximumNs == null) {
      return true;
    }
    final int size = s.length();
    int ncountLeft = 0;
    int ncountRight = 0;
    final int half = size / 2;
    if (half > mMaxLength) {
      mMaxLength = half;
    }
    if (half < mMinLength) {
      mMinLength = half;
    }
    for (int k = 0, j = half; k < half; k++, j++) {
      if (Character.toLowerCase(s.charAt(k)) == 'n') {
        ncountLeft++;
      }
      if (Character.toLowerCase(s.charAt(j)) == 'n') {
        ncountRight++;
      }
    }
    final boolean r = ncountLeft <= mMaximumNs && ncountRight <= mMaximumNs;
    if (r) {
      if (ncountLeft > mMaximumNsCount) {
        mMaximumNsCount = ncountLeft;
      }
      if (ncountRight > mMaximumNsCount) {
        mMaximumNsCount = ncountRight;
      }
    } else {
      mSkippedNReads++;
      mSkippedNResidues += size;
    }
    return r;
  }

  private boolean processLine() throws IOException {
    String line;
    if (mReader == null) {
      mStream = new FileInputStream(mReads);
      final InputStream gi;
      if (FileUtils.isGzipFilename(mReads)) {
        try {
          gi = new AsynchInputStream(GzipUtils.createGzipInputStream(mStream)); // BufferedInputStream doesn't seem to help here
        } catch (final IOException e) {
          throw new NoTalkbackSlimException(ErrorType.NOT_A_CG_INPUT, "File not in GZIP format");
        }
      } else if (FileUtils.isBzip2Filename(mReads)) {
         try {
           gi = new AsynchInputStream(new CBZip2InputStream(new BufferedInputStream(mStream, FileUtils.BUFFERED_STREAM_SIZE)));
        } catch (final IOException e) {
          throw new NoTalkbackSlimException(ErrorType.NOT_A_CG_INPUT, "File not in BZIP2 format");
        }
      } else {
        gi = mStream;
      }
      mReader = new BufferedReader(new InputStreamReader(gi));
    }
    while ((line = mReader.readLine()) != null) {
      if (line.length() > 0 && line.charAt(0) != '>' && line.charAt(0) != '#') {
        final String[] parts = PATTERN_TAB_SPLIT.split(line);
        if (parts.length != 3) {
          Diagnostic.userLog("Malformed input for sequence " + mSeqNo + ": " + line);
          throw new NoTalkbackSlimException(ErrorType.NOT_A_CG_INPUT, "File contains invalid CG data. Insufficient number of fields.");
        }
        if (parts[READ_FIELD].length() != parts[QUALITY_FIELD].length()) {
          Diagnostic.userLog("Malformed input for sequence " + mSeqNo + ": " + line);
          throw new NoTalkbackSlimException(ErrorType.NOT_A_CG_INPUT, "File contains invalid CG data. Mismatched base/quality field lengths.");
        }
        // parts[0] is flags, ignore for now
        if (acceptableNs(parts[READ_FIELD])) {
          mCurrentLength = parts[READ_FIELD].length() / 2;

          if (mLeftRead == null || mLeftRead.length < mCurrentLength) {
            mLeftRead = new byte[mCurrentLength];
            mRightRead = new byte[mCurrentLength];
            mLeftQuality = new byte[mCurrentLength];
            mRightQuality = new byte[mCurrentLength];
          }

          final byte[] readBytes = parts[READ_FIELD].getBytes();
          copySequenceDataField(readBytes, 0, mLeftRead, 0, mCurrentLength);
          copySequenceDataField(readBytes, mCurrentLength, mRightRead, 0, mCurrentLength);
          copyQualityDataField(parts[QUALITY_FIELD], 0, mLeftQuality, 0, mCurrentLength);
          copyQualityDataField(parts[QUALITY_FIELD], mCurrentLength, mRightQuality, 0, mCurrentLength);

          mSeqNo++;
          return true;
        }
      }
    }
    return false;
  }

  private void copySequenceDataField(byte[] src, int srcPos, byte[] dest, int destPos, int length) {
    byte residueByte;
    for (int i = srcPos, j = destPos; i < srcPos + length; i++, j++) {
      residueByte = mFastaSymbolLookupTable[src[i]];
      if (residueByte == (byte) 255) {
        //unrecognized character, print warning, shove in unknown
        Diagnostic.warning(WarningType.BAD_TIDE, name(), Character.toString((char) src[i]), CHAR_TO_RESIDUE.unknownResidue().toString());
        residueByte = (byte) CHAR_TO_RESIDUE.unknownResidue().ordinal();
      }
      dest[j] = residueByte; //(byte) tempResidue.ordinal();
    }
  }

  private void copyQualityDataField(String src, int srcPos, byte[] dest, int destPos, int length) {
    char tempChar;
    for (int i = srcPos, j = destPos; i < srcPos + length; i++, j++) {
      tempChar = src.charAt(i);
      dest[j] = (byte) (tempChar - '!');
    }
  }

  // First call to nextSequence() will make this true
  private boolean mLeft = false;

  @Override
  public String name() {
    // Original file has no names, so we just make some up
    return mSeqNo + "-" + (mLeft ? "A" : "B");
  }

  @Override
  public long getWarningCount() {
    return 0;
  }

  @Override
  public boolean nextSequence() throws IOException {
    mLeft = !mLeft;
    if (mLeft) {
      final boolean ok = processLine();
      if (!ok) {
        return false;
      }
    }
    return true;
  }

  private static final DNAFastaSymbolTable CHAR_TO_RESIDUE = new DNAFastaSymbolTable();

  @Override
  public byte[] sequenceData() {
    return mLeft ? mLeftRead : mRightRead;
  }

  /**
   * Return Phred scores.
   * @return integers between 0 and 63
   */
  @Override
  public byte[] qualityData() {
    return mLeft ? mLeftQuality : mRightQuality;
  }

  /**
   * Not implemented
   */
  @Override
  public void setDusting(final boolean val) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int getMaxNCount() {
    return mMaximumNsCount;
  }

  @Override
  public int getSkippedReads() {
    return mSkippedNReads;
  }

  @Override
  public long getSkippedResidues() {
    return mSkippedNResidues;
  }

  @Override
  public long getMaxLength() {
    return mMaxLength;
  }

  @Override
  public long getMinLength() {
    return mMinLength;
  }

  @Override
  public int currentLength() {
    return mCurrentLength;
  }

  @Override
  public long getDusted() {
    return 0;
  }
}
