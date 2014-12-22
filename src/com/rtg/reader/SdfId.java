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

import java.util.UUID;

import com.rtg.util.StringUtils;

/**
 * Encapsulate <code>SDF-ID</code>'s in its various incarnations
 */
public class SdfId {

  private final long mSdfId;
  private final UUID mSdfUuid;
  private final boolean mHasUuid;

  private final int mLongHash;

  /**
   * Create a new random UUID style SDF-ID
   */
  public SdfId() {
    this(UUID.randomUUID());
  }

  /**
   * Create a UUID style SDF-ID
   * @param id the id
   */
  public SdfId(UUID id) {
    mSdfUuid = id;
    mSdfId = 0;
    mHasUuid = true;
    mLongHash = Long.valueOf(mSdfId).hashCode();
  }

  /**
   * Create an old long based SDF-ID
   * @param id the id
   */
  public SdfId(long id) {
    mSdfUuid = new UUID(0, 0);
    mSdfId = id;
    mHasUuid = false;
    mLongHash = Long.valueOf(mSdfId).hashCode();
  }

  /**
   * <p>Construct an SDF-ID from given string, the type is determined by the String itself
   * <p> i.e. <code>xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx</code> gives a UUID style SDF-ID
   * <p> i.e. <code>xxxxxxxxxxxxxxxx</code> gives a long based SDF-ID
   * <p><code>x</code> means a hexadecimal digit
   * @param value value to decode
   */
  public SdfId(String value) {
    if (value.contains("-")) {
      mHasUuid = true;
      mSdfUuid = UUID.fromString(value);
      mSdfId = 0;
    } else {
      mHasUuid = false;
      mSdfUuid = new UUID(0, 0);
      mSdfId = StringUtils.fromLongHexString(value);
    }
    mLongHash = Long.valueOf(mSdfId).hashCode();
  }

  /**
   * Checks if this <code>SDF-ID</code> is the same as the supplied one. Except when
   * one or the other has no id (i.e. long constructor with 0 as argument), this implies
   * that no id was available and we cannot check for equality, either because the SDF has
   * no id or that version of the output files did not record the id.
   * @param other other id to check
   * @return true if both <code>SdfId</code>'s are the same or one is not available. false otherwise.
   */
  public boolean check(SdfId other) {
    if (mHasUuid) {
      if (other.mHasUuid) {
        return mSdfUuid.equals(other.mSdfUuid);
      } else {
        return other.mSdfId == 0;
      }
    }
    //mHasUuid == false
    if (other.mHasUuid) {
      return mSdfId == 0;
    }
    //mHasUuid == false && other.mHasUuid == false
    return mSdfId == other.mSdfId || mSdfId == 0 || other.mSdfId == 0;
  }

  @Override
  public String toString() {
    if (mHasUuid) {
      return mSdfUuid.toString();
    } else {
      return Long.toHexString(mSdfId);
    }
  }

  long getLowBits() {
    if (mHasUuid) {
      return mSdfUuid.getLeastSignificantBits();
    }
    return mSdfId;
  }

  long getHighBits() {
    if (mHasUuid) {
      return mSdfUuid.getMostSignificantBits();
    }
    //for testing only
    return 0L;
    //throw new UnsupportedOperationException("High bits only available from UUID type SDF ID");
  }

  /**
   * Reports if a valid SDF-ID in contained within, non-valid SDF-IDs are constructed
   * by calling the {@link SdfId#SdfId(long)} constructor with 0 as the argument. This is used to indicate
   * an SDF-ID was not available from whatever source.
   * @return true if this is a valid SDF-ID
   */
  public boolean available() {
    return mHasUuid || mSdfId != 0;
  }

  @Override
  public boolean equals(Object obj) {
    if (obj instanceof SdfId) {
      final SdfId other = (SdfId) obj;
      if (mHasUuid != other.mHasUuid) {
        return false;
      }
      if (mHasUuid) {
        return mSdfUuid.equals(other.mSdfUuid);
      }
      return mSdfId == other.mSdfId;
    }
    return false;
  }

  @Override
  public int hashCode() {
    if (mHasUuid) {
      return mSdfUuid.hashCode();
    }
    return mLongHash;
  }
}
