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
public final class SimulatedReadNameParserFactory {

  /** Private to prevent instantiation. */
  private SimulatedReadNameParserFactory() { }

  /**
   * Chooses a suitable parser implementation for the supplied simulated read name.
   *
   * @param exampleReadName an example read name
   * @return the parser, or null if the supplied name does not appear to be a handled format.
   */
  public static SimulatedReadNameParser getParser(String exampleReadName) {
    final SimulatedReadNameParser parser;
    if (NewestReadNameParser.looksOk(exampleReadName)) {
      parser = new NewestReadNameParser();
    } else if (NewReadNameParser.looksOk(exampleReadName)) {
      parser = new NewReadNameParser();
    } else if (OldReadNameParser.looksOk(exampleReadName)) {
      parser = new OldReadNameParser();
    } else if (DwgsimReadNameParser.looksOk(exampleReadName)) {
      parser = new DwgsimReadNameParser();
    } else {
      parser = null;
    }
    if (parser != null) {
      try {
        if (parser.setReadInfo(exampleReadName, 0)) {
          return parser;
        }
      } catch (RuntimeException e) {
        return null;
      }
    }
    return null;
  }

  /**
   * Chooses a suitable parser implementation for the supplied simulated read name.
   * @param exampleReadName an example read name
   * @param offsets the offsets for the samples in the mutations VCF file to apply corrections for.
   * @return the parsers, or null if the supplied name does not appear to be a handled format.
   */
  public static SimulatedReadNameParser[] getParsers(String exampleReadName, MutatedSampleOffsets... offsets) {
    final SimulatedReadNameParser[] parsers = new SimulatedReadNameParser[offsets.length];
    for (int i = 0; i < exampleReadName.length(); ++i) {
      for (int j = 0; j < offsets.length; ++j) {
        if (offsets[j] != null) {
          final SimulatedReadNameParser p = getParser(exampleReadName);
          if (p != null) {
            parsers[j] = new MutatedReferenceReadNameParser(p, offsets[j]);
          }
        }
      }
    }
    return parsers;
  }

}
