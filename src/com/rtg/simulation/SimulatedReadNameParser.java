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
package com.rtg.simulation;

/**
 * Interface to represent parsers of simulated read names
 */
public interface SimulatedReadNameParser {

  /**
   * Sets the read name information to parse.
   *
   * @param readName the read name.
   * @param readLength the read length.
   * @return false if the read name does not appear to be a simulated read name that can be parsed
   */
  boolean setReadInfo(String readName, int readLength);

  /**
   * @return true if the read appears to be generated as a chimera
   */
  boolean isChimera();

  /**
   * @return true if the read appears to be generated as a PCR or optical duplicate
   */
  boolean isDuplicate();

  /**
   * Return the one based position on the template from which the read was generated.
   * @return the position
   */
  long templatePosition();

  /**
   * Return the name of the template sequence from which the read was generated.
   * @return the name
   */
  String templateName();

  /**
   * Return the id of the template set (i.e. which haploid set, when
   * doing diploid simulation) from which the read was generated.
   * @return the id of the template set
   */
  int templateSet();

  /**
   * Returns true if the read was generated in the forward frame
   * @return true if forward frame
   */
  boolean forwardFrame();

  /**
   * @return the unique id of the read
   */
  int readId();

  /**
   * Return the unique label of the read, with any meta-data stripped
   * @return the name
   */
  String readName();

  /**
   * Return the number of substitutions made in this read.
   * @return the number of substitutions
   */
  int substitutions();

  /**
   * Return the number of insertions made in this read.
   * @return the number of insertions
   */
  int insertions();

  /**
   * Return the number of deletions made in this read.
   * @return the number of deletions
   */
  int deletions();

  /**
   * @return the total number of mismatches made in this read (comparable to <code>NM</code> attribute in SAM record)
   */
  int numMismatches();
}
