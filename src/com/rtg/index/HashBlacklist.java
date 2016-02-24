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
package com.rtg.index;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Load a blacklist
 */
public final class HashBlacklist {


  private HashBlacklist() { }

  private static final String BLACKLIST_SUBDIR = "blacklists";

  /**
   * Load a blacklist file (assumed to be lines of "k-mer[TAB]count")
   * @param sdfDir directory containing sdf
   * @return list of hashes
   * @throws IOException if an IO error occurs
   */
  public static List<Long> loadBlacklist(File sdfDir, int wordSize, int threshold) throws IOException {
    final ArrayList<Long> ret = new ArrayList<>();
    final File blacklistFile = new File(sdfDir, BLACKLIST_SUBDIR + File.separator + "w" + wordSize);
    Diagnostic.developerLog("Loading blacklist at word size " + wordSize + " and threshold " + threshold);
    try (BufferedReader br = new BufferedReader(new FileReader(blacklistFile))) {
      String line;
      while ((line = br.readLine()) != null) {
        final String[] split = line.split("\t");
        final String word = split[0];
        final long count = Long.parseLong(split[1]);
        if (count >= threshold) {
          final byte[] dna = DnaUtils.encodeArray(word.getBytes());
          ret.add(dnaToLong(dna));
        }
      }
    }
    Diagnostic.developerLog("Loaded " + ret.size() + " hashes from blacklist");
    return ret;
  }

  private static long dnaToLong(byte[] dna) {
    long ret = 0;
    for (int i = 0; i < dna.length; i++) {
      ret <<= 2;
      ret |= (dna[i] - 1) & 0b11L;
    }
    return ret;
  }
}
