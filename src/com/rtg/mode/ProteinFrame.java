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
package com.rtg.mode;

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

/**
 * Protein can only have a single untranslated frame (no translation possible
 * and it makes no sense to take a reverse complement).
 */
public final class ProteinFrame implements Frame, PseudoEnum {
  /** Protein. */
  public static final ProteinFrame PROTEIN = new ProteinFrame(0, "PROTEIN");

  private final int mOrdinal;
  private final String  mName;

  private ProteinFrame(final int ordinal, final String name) {
    mOrdinal = ordinal;
    mName = name;
  }

  @Override
  public String toString() {
    return mName;
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public int ordinal() {
    return mOrdinal;
  }

  private static final EnumHelper<ProteinFrame> HELPER = new EnumHelper<>(ProteinFrame.class, new ProteinFrame[] {PROTEIN});

  /**
   * @return the list of enum values
   */
  public static ProteinFrame[] values() {
    return HELPER.values();
  }

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str name of enum
   * @return the enum value
   */
  public static ProteinFrame valueOf(final String str) {
    return HELPER.valueOf(str);
  }

  @Override
  public String display() {
    return "";
  }

  @Override
  public boolean isForward() {
    return true;
  }

  @Override
  public Frame getReverse() {
    throw new UnsupportedOperationException("Not supported");
  }

  /**
   * Get the frame corresponding to the integer value.
   * @param value int value
   * @return the frame
   */
  static ProteinFrame frameFromCode(final int value) {
    if (value != 0) {
      throw new IllegalArgumentException(String.valueOf(value));
    }
    return PROTEIN;
  }

  private static final byte UNKNOWN_RESIDUE = (byte) Protein.X.ordinal();

  @Override
  public byte code(final byte[] codes, final int length, final int index, int offset, int fullLength) {
    return index >= length ? UNKNOWN_RESIDUE : codes[index];
  }

  @Override
  public byte code(byte[] codes, int length, int index) {
    return code(codes, length, index, 0, length);
  }
  @Override
  public int phase() {
    return 0;
  }

  @Override
  public int calculateFirstValid(int offset, int length, int fullLength) {
    return 0;
  }

  @Override
  public int calculateLastValid(int offset, int length, int fullLength) {
    return length;
  }

}

