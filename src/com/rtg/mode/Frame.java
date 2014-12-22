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

/**
 * The possible framing for sequences.
 *
 */
public interface Frame {

  /**
   * Get a unique integer code for the frame.
   * @return a unique integer code.
   */
  int ordinal();

  /**
   * Get a code with respect to the frame.
   * @param codes the underlying array of codes.
   * @param length the number of valid entries in codes.
   * @param index the position to access with respect to this frame.
   *              note that it is in the units of codes not of the translated
   *              results - for example in translated frames index will typically
   *              be incremented by 3 to get succesive translated codes.
   * @return the code at the position specified by index.
   */
  byte code(byte[] codes, int length, int index);

  /**
   * Get a code with respect to the frame.
   * @param codes the underlying array of codes.
   * @param length the number of valid entries in codes.
   * @param index the position to access with respect to this frame.
   *              note that it is in the units of codes not of the translated
   *              results - for example in translated frames index will typically
   *              be incremented by 3 to get succesive translated codes.
   * @param firstValid start position of first valid code in code array material (inclusive).
   * @param lastValid last valid code in code array material (exclusive). for a translated frame this will be the
   * position of the first code that cannot be included in a valid translated result.
   * @return the code at the position specified by index.
   */
  byte code(byte[] codes, int length, int index, int firstValid, int lastValid);

  /**
   * Calculate the first valid code. (inclusive)
   * @param offset position into source material that array starts
   * @param length length of data in array
   * @param fullLength length of source material
   * @return the first valid code to pass to code method
   */
  int calculateFirstValid(int offset, int length, int fullLength);

  /**
   * Calculate the last valid code. (exclusive)
   * @param offset position into source material that array starts
   * @param length length of data in array
   * @param fullLength length of source material
   * @return the last valid code to pass to code method
   */
  int calculateLastValid(int offset, int length, int fullLength);


  /**
   * Get the string to be used in output formats for describing the frame.
   * @return a string to be used in output formats.
   */
  String display();

  /**
   * Get the offset of the frame with respect to the underlying sequence.
   * Is &gt;= 0 (0 except for translated frames).
   * @return the offset of the frame with respect to the underlying sequence.
   */
  int phase();

  /**
   * Check if the frame is forward (or untranslated).
   * @return true iff frame is forward.
   */
  boolean isForward();

  /**
   * Get the reverse of this frame
   * @return a frame
   */
  Frame getReverse();
}

