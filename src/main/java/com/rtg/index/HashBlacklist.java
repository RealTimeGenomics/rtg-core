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
    Diagnostic.userLog("Loading blacklist at word size " + wordSize + " and threshold " + threshold);
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
    for (final byte base : dna) {
      ret <<= 2;
      ret |= (base - 1) & 0b11L;
    }
    return ret;
  }
}
