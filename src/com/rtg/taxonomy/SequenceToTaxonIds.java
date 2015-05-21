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
package com.rtg.taxonomy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;

import com.rtg.util.StringUtils;

/**
 */
public final class SequenceToTaxonIds {

  private SequenceToTaxonIds() { }

  /**
   * Load a map from (short) sequence names to taxon ids supplied as tab separated input
   * @param tsv the input file
   * @return the map
   * @throws IOException if there are problems reading the input
   */
  public static Map<String, Integer> sequenceToIds(File tsv) throws IOException {
    return sequenceToIds(new FileReader(tsv));
  }

  /**
   * Load a map from (short) sequence names to taxon ids supplied as tab separated input
   * @param input input source
   * @return the map
   * @throws IOException if there are problems reading the input
   */
  public static Map<String, Integer> sequenceToIds(Reader input) throws IOException {
    final Map<String, Integer> ret = new HashMap<>();
    try (BufferedReader mergedReader = new BufferedReader(input)) {
      String line;
      while ((line = mergedReader.readLine()) != null) {
        if (line.startsWith("#") || line.length() == 0) {
          continue;
        }
        final String[] parts = StringUtils.split(line, '\t');
        if (parts.length != 2) {
          throw new IOException("Malformed line: " + line);
        }
        final String name;
        final int indexOfSpace = parts[1].indexOf(" ");
        if (indexOfSpace != -1) {
          name = parts[1].substring(0, indexOfSpace);
        } else {
          name = parts[1];
        }

        final int taxId;
        try {
          taxId = Integer.parseInt(parts[0]);
        } catch (NumberFormatException e) {
          throw new IOException("Malformed line: " + line);
        }
        if (ret.put(new String(name.toCharArray()), taxId) != null) {
          throw new IOException("Duplicate name detected: " + line) ;
        }
      }
    }
    return ret;
  }
}
