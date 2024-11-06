/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
