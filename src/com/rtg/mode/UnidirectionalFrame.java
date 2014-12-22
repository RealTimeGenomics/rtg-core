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
 * A nucleotide sequence that is not translated.
 */
public final class UnidirectionalFrame implements Frame, PseudoEnum {

  /**
   * Normal forward direction.
   */
  public static final UnidirectionalFrame FORWARD = new UnidirectionalFrame(0, "FORWARD");

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

  private int mOrdinal;
  private String mName;

  private UnidirectionalFrame(final int ordinal, final String name) {
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

  private static final EnumHelper<UnidirectionalFrame> HELPER = new EnumHelper<>(UnidirectionalFrame.class, new UnidirectionalFrame[] {FORWARD});

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str name of value
   * @return the enum value
   */
  public static UnidirectionalFrame valueOf(final String str) {
    return HELPER.valueOf(str);
  }

  /**
   * @return list of enum values
   */
  public static UnidirectionalFrame[] values() {
    return HELPER.values();
  }

  /**
   * Get the frame corresponding to the integer value.
   * @param value int value
   * @return the frame
   */
  static UnidirectionalFrame frameFromCode(final int value) {
    return values()[value];
  }

  @Override
  public byte code(final byte[] codes, final int length, final int index, int offset, int fullLength) {
    if (index >= length) {
      throw new RuntimeException("length=" + length + " index=" + index);
    }
    return codes[index];
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

