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
package com.rtg.ngs.tempstage;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.ReadableByteChannel;

import com.rtg.sam.SamBamConstants;

/**
 * A Binary writable representation of a temp file record
 */
public class BinaryTempFileRecord {

  private static final int PAIRED_FIELDS = 0x1;
  private static final int MD_FIELDS = 0x2;
  private static final int CG_FIELDS = 0x4;
  private static final int UNFILTERED_MATE_FIELDS = 0x8;

  private final byte mFields;

  /**
   * Create a binary temp file record
   * @param paired true if this is paired end data
   * @param legacyCigars true if legacy cigar mode is enabled
   * @param cg true if the data is from Complete Genomics
   * @param unfiltered true if all hits mode is enabled
   */
  public BinaryTempFileRecord(boolean paired, boolean legacyCigars, boolean cg, boolean unfiltered) {
    int f = 0;
    f |= paired ? PAIRED_FIELDS : 0;
    f |= legacyCigars ? MD_FIELDS : 0;
    f |= cg ? CG_FIELDS : 0;
    f |= unfiltered && paired ? UNFILTERED_MATE_FIELDS : 0;
    mFields = (byte) f;
  }

  //Fixed fields
  int mReferenceId;
  int mReadId;
  byte mSamFlags;
  int mStartPosition;
  byte[] mCigarString = null;
  int mAlignmentScore;
  int mNumberMismatches;

  //MD field
  byte[] mMdString = null;

  //Paired end fields
  int mMatePosition;
  int mTemplateLength;
  int mComboScore;
  //CG fields
  byte[] mReadString;
  byte[] mSuperCigarString;
  byte[] mReadDeltaString;

  //unfiltered mating fields
  private boolean mUnfilteredMated;

  /**
   * @return True iff this is the phony record to mark no more records
   */
  public boolean isSentinelRecord() {
    return mReferenceId == -1;
  }

  /**
   * Mark this record as end of file marker. No other values should be set.
   */
  public void setSentinelRecord() {
    mReferenceId = -1;
    mCigarString = new byte[0];
    if (hasMdField()) {
      mMdString = new byte[0];
    }
    if (hasCgFields()) {
      mReadString = new byte[0];
      mSuperCigarString = new byte[0];
      mReadDeltaString = new byte[0];
    }
  }

  public void setReferenceId(int referenceId) {
    mReferenceId = referenceId;
  }

  /**
   * @param readId internal <code>readId</code>, should be between 0 and Integer.MAX_VALUE
   */
  public void setReadId(int readId) {
    mReadId = readId;
  }

  public void setSamFlags(byte samFlags) {
    mSamFlags = samFlags;
  }

  /**
   * @param startPosition one based alignment start position
   */
  public void setStartPosition(int startPosition) {
    mStartPosition = startPosition;
  }

  public void setCigarString(byte[] cigarString) {
    mCigarString = cigarString;
  }

  public void setAlignmentScore(int alignmentScore) {
    mAlignmentScore = alignmentScore;
  }

  public void setNumberMismatches(int numberMismatches) {
    mNumberMismatches = numberMismatches;
  }

  public void setMdString(byte[] mdString) {
    mMdString = mdString;
  }

  public int getAlignmentScore() {
    return mAlignmentScore;
  }

  public byte[] getCigarString() {
    return mCigarString;
  }

  public byte[] getMdString() {
    return mMdString;
  }


  public int getNumberMismatches() {
    return mNumberMismatches;
  }


  /**
   * @return internal <code>readId</code>, should be between 0 and Integer.MAX_VALUE
   */
  public int getReadId() {
    return mReadId;
  }


  public int getReferenceId() {
    return mReferenceId;
  }


  public byte getSamFlags() {
    return mSamFlags;
  }


  /**
   * @return one based alignment start position
   */
  public int getStartPosition() {
    return mStartPosition;
  }

  public boolean isReverseStrand() {
    return (mSamFlags & SamBamConstants.SAM_READ_IS_REVERSE) != 0;
  }

  /**
   * @return true if this read is paired during sequencing
   */
  public boolean isReadPaired() {
    return (mSamFlags & SamBamConstants.SAM_READ_IS_PAIRED) != 0;
  }

    /**
   * @param matePosition one based start position of mates alignment
   */
  public void setMatePosition(int matePosition) {
    mMatePosition = matePosition;
  }

//  public void setMateReferenceId(int mateReferenceId) {
//    mMateReferenceId = mateReferenceId;
//  }

  public void setTemplateLength(int templateLength) {
    mTemplateLength = templateLength;
  }

  public void setComboScore(int comboScore) {
    mComboScore = comboScore;
  }

  /**
   * @return one based start position of mates alignment
   */
  public int getMatePosition() {
    return mMatePosition;
  }

//  public int getMateReferenceId() {
//    return mMateReferenceId;
//  }

  /**
   * @return the length this pair on the template (from read start of leftmost read on template, to read end of rightmost read on template)
   */
  public int getTemplateLength() {
    return mTemplateLength;
  }

  public int getComboScore() {
    return mComboScore;
  }

  public boolean isMateReverseStrand() {
    return (mSamFlags & SamBamConstants.SAM_MATE_IS_REVERSE) != 0;
  }

  public boolean isFirstOfPair() {
    return (mSamFlags & SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR) != 0;
  }

  public boolean isSecondOfPair() {
    return (mSamFlags & SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR) != 0;
  }

  public byte[] getCgReadString() {
    return mReadString;
  }

  public void setCgReadString(byte[] readString) {
    mReadString = readString;
  }

  public byte[] getSuperCigarString() {
    return mSuperCigarString;
  }

  public void setSuperCigarString(byte[] superCigarString) {
    mSuperCigarString = superCigarString;
  }

  /**
   * @return The read delta (a.k.a. XR field in sam)
   */
  public byte[] getReadDeltaString() {
    return mReadDeltaString;
  }

  public void setReadDeltaString(byte[] readDeltaString) {
    mReadDeltaString = readDeltaString;
  }


  public boolean isUnfilteredMated() {
    return mUnfilteredMated;
  }

  public void setUnfilteredMated(boolean mated) {
    mUnfilteredMated = mated;
  }

  /**
   * Write this alignment record to a data output stream.
   * @param out the data output stream
   */
  public void writeNio(ByteBuffer out) {
    if (isSentinelRecord()) {
      out.putInt(-1);
      return;
    }
    out.putInt(mReferenceId);
    out.putInt(mStartPosition);
    out.putInt(mReadId);
    out.put(mSamFlags);
    out.putInt(mCigarString.length);
    out.put(mCigarString);
    out.putInt(mAlignmentScore);
    out.putInt(mNumberMismatches);
    if (hasMdField()) {
      out.putInt(mMdString.length);
      out.put(mMdString);
    }

    if (hasPairedFields()) {
      out.putInt(mMatePosition);
      out.putInt(mTemplateLength);
      out.putInt(mComboScore);
    }
    if (hasCgFields()) {
      writeArrayNio(out, mReadString);
      writeArrayNio(out, mSuperCigarString);
      writeArrayNio(out, mReadDeltaString);
    }
    if (hasUnfilteredField()) {
      out.put((byte) (mUnfilteredMated ? 1 : 0));
    }
  }

  /**
   * Read into this alignment record from a data input stream.
   * @param in the data input stream
   * @param ch channel to read further data from if necessary
   * @throws IOException if an exception occurs while reading
   * @return true if a record was read, false if the end of the input data has been reached
   */
  public boolean readNio(ByteBuffer in, ReadableByteChannel ch) throws IOException {
    if (in.remaining() < 100) {
      readMore(in, ch);
    }
    mReferenceId = in.getInt();
    if (mReferenceId == -1) {
      return false;
    }
    mStartPosition = in.getInt();
    mReadId = in.getInt();
    mSamFlags = in.get();
    final int cigarLen = in.getInt();
    if (in.remaining() < cigarLen + 12) {
      readMore(in, ch);
    }
    mCigarString = new byte[cigarLen];
    in.get(mCigarString);
    mAlignmentScore = in.getInt();
    mNumberMismatches = in.getInt();
    if (hasMdField()) {
      final int mdStringLen = in.getInt();
      mMdString = new byte[mdStringLen];
      if (in.remaining() < mdStringLen) {
        readMore(in, ch);
      }
      in.get(mMdString);
    }

    if (hasPairedFields()) {
      if (in.remaining() < 12) {
        readMore(in, ch);
      }
      mMatePosition = in.getInt();
      mTemplateLength = in.getInt();
      mComboScore = in.getInt();
    }

    if (hasCgFields()) {
      mReadString = readArrayNio(in, ch);
      mSuperCigarString = readArrayNio(in, ch);
      mReadDeltaString = readArrayNio(in, ch);
    }
    if (hasUnfilteredField()) {
      if (in.remaining() < 1) {
        readMore(in, ch);
      }
      final byte t = in.get();
      mUnfilteredMated = t != 0;
    }
    return true;
  }

  static void readMore(ByteBuffer in, ReadableByteChannel ch) throws IOException {
    in.compact();
    int rr = 0;
    while (in.remaining() != 0 && rr >= 0) {
      rr = ch.read(in);
    }
    in.flip();
  }

  private static void writeArrayNio(ByteBuffer out, byte[] arr) {
    out.putInt(arr.length);
    out.put(arr);
  }
  private static byte[] readArrayNio(ByteBuffer in, ReadableByteChannel ch) throws IOException {
    if (in.remaining() < 4) {
      readMore(in, ch);
    }
    final int length = in.getInt();
    if (in.remaining() < length) {
      readMore(in, ch);
    }
    final byte[] ret = new byte[length];
    in.get(ret);
    return ret;
  }

  private boolean hasPairedFields() {
    return (mFields & PAIRED_FIELDS) != 0;
  }
  private boolean hasMdField() {
    return (mFields & MD_FIELDS) != 0;
  }
  private boolean hasCgFields() {
    return (mFields & CG_FIELDS) != 0;
  }
  private boolean hasUnfilteredField() {
    return (mFields & UNFILTERED_MATE_FIELDS) != 0;
  }
}
