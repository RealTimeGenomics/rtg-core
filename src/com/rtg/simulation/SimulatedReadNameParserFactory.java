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
    for (int i = 0; i < exampleReadName.length(); i++) {
      for (int j = 0; j < offsets.length; j++) {
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
