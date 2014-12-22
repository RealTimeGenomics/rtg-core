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
 * Represents one monomer in a biopolymer.
 * Likely implementations include nucleotides and polypeptides.
 *
 */
public interface Residue {

  /**
   * Get a numeric code for the object.
   * It may be assumed that the codes range contiguously from 0 to a maximum
   * value and that any unknown residue accurs at the end of this sequence and that
   * there is precisely one unknown residue.
   * @return numeric code for object.
   */
  int ordinal();

  /**
   * Check if this residue is to be ignored when forming a window.
   * @return true if this residue is to be ignored.
   */
  boolean ignore();

  /**
   * Get the type of sequence this residue will be part of.
   * @return the sequence type
   */
  SequenceType type();

  /**
   * Get the single letter code for this residue.
   * @return the single letter code for this residue.
   */
  String toString();
}

