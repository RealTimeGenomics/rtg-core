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
package com.rtg.variant.realign;


/**
 * Stores information about a read, a template and their relationship.
 *
 */
public interface Environment {

  /**
   * Get the probability of the specified nucleotide
   * at the given index in the read. These probabilities will typically
   * be computed from quality scores.
   * @param index position in read (0 based).
   * @return the raw probability (not its log).
   */
  double quality(int index);

  /**
   * Get the base at specified position in the read.
   * @param index position in read (0 based).
   * @return the nucleotide (0 = N).
   */
  byte read(int index);

  /**
   * Get the base at the specified relative position in the template.
   * @param index position in template (0 based from the expected start of the read).
   * @return the nucleotide (0 = N).
   */
  byte template(int index);

  /**
   * Convert a relative position on the template to its real position.
   * This does NOT take into account the actual alignment of the read,
   * or any CG gaps.  It simply assumes that the read is laid along the
   * template one-for-one.
   * @param index zero means the nominal 'start' of the read on the template
   *   (which will be the end of the read for a reversed read).
   * @return the real template position (may be off the template).
   */
  int absoluteTemplatePosition(int index);

  /**
   * Get a position which is guaranteed to be greater than the position of the last non-null nucleotide in the read
   * (0 based).
   * @return the length of the read.
   */
  int readLength();

  /**
   * @return the length of the template.
   */
  int templateLength();

  /**
   * Get the maximum shift permitted when doing alignment.
   * @return the maximum shift permitted when doing alignment.
   */
  int maxShift();
}
