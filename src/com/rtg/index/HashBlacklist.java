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
import java.io.Reader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Load a blacklist
 */
public final class HashBlacklist {


  private HashBlacklist() { }

  static final String BLACKLIST_SUBDIR = "blacklists";

  /**
   * Determine if SDF has a blacklist for the appropriate word size
   * @param sdfDir directory containing SDF
   * @param wordSize kmer size
   * @return true if blacklist exists
   */
  public static boolean blacklistExists(File sdfDir, int wordSize) {
    final File blacklist = getFile(sdfDir, wordSize);
    return blacklist.exists();
  }

  /**
   * Load a blacklist file (assumed to be lines of "k-mer[TAB]count")
   * @param sdfDir directory containing sdf
   * @param wordSize kmer size
   * @param threshold only load hashes with count meeting this threshold
   * @return list of hashes
   * @throws IOException if an IO error occurs
   */
  public static List<Long> loadBlacklist(File sdfDir, int wordSize, int threshold) throws IOException {
    final File blacklistFile = getFile(sdfDir, wordSize);
    Diagnostic.developerLog("Loading blacklist at word size " + wordSize + " and threshold " + threshold);
    return loadBlacklist(new FileReader(blacklistFile), threshold);
  }

  /**
   * Installs given blacklist into given SDF
   * @param blacklist blacklist to install
   * @param sdfDir SDF to install into
   * @param wordSize kmer size
   * @throws IOException If blacklist already exists or other IO error
   */
  public static void installBlacklist(File blacklist, File sdfDir, int wordSize) throws IOException {
    final File destination = getFile(sdfDir, wordSize);
    if (destination.exists()) {
      throw new IOException("Blacklist already exists in " + sdfDir + " for word size " + wordSize);
    }
    if (!destination.getParentFile().isDirectory() && !destination.getParentFile().mkdir()) {
      throw new IOException("Could not make blacklist directory");
    }
    Files.copy(blacklist.toPath(), destination.toPath());
  }

  /**
   * load blacklist from reader
   * @param blacklistData blacklist data (assumed to be lines of "k-mer[TAB]count")
   * @param threshold only load hashes with count meeting this threshold
   * @return list of hashes
   * @throws IOException if an IO error occurs
   */
  public static List<Long> loadBlacklist(Reader blacklistData, int threshold) throws IOException {
    final ArrayList<Long> ret = new ArrayList<>();
     try (BufferedReader br = new BufferedReader(blacklistData)) {
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

  private static File getFile(File sdfDir, int wordSize) {
    return new File(sdfDir, BLACKLIST_SUBDIR + File.separator + "w" + wordSize);
  }

  /**
   * @param wordSize kmer size
   * @return number of bits used in hash for kmer
   */
  public static int hashBits(int wordSize) {
    return 2 * wordSize;
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
