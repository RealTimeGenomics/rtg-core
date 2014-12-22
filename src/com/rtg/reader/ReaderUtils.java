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
package com.rtg.reader;


import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Utility functions for manipulating files that are not provided in the File
 * class.
 *
 */
public final class ReaderUtils {

  private ReaderUtils() {
  }

  /**
   * Tells you if the directory passed in contains a paired end SDF left / right folder
   *
   * @param dir the input directory
   * @return true if the left / right folders exist within the given directory; otherwise false
   */
  public static boolean isPairedEndDirectory(final File dir) {
    if (dir == null || !dir.exists() || !dir.isDirectory()) {
      return false;
    }
    final File lDir = new File(dir, "left");
    final File rDir = new File(dir, "right");
    return lDir.exists() && rDir.exists() && lDir.isDirectory() && rDir.isDirectory();
  }

  /**
   * Gets the left end of an existing paired end SDF
   *
   * @param dir the input directory
   * @return the left end SDF directory
   * @throws CorruptSdfException if one of the directories has been deleted since the check
   */
  public static File getLeftEnd(final File dir) throws CorruptSdfException {
    if (dir == null || !dir.exists() || !dir.isDirectory()) {
      throw new CorruptSdfException(dir);
    }
    final File lDir = new File(dir, "left");
    if (lDir.exists() && lDir.isDirectory()) {
      return lDir;
    }
    throw new CorruptSdfException(dir);
  }

  /**
   * Gets the right end of an existing paired end SDF
   *
   * @param dir the input directory
   * @return the right end SDF directory
   * @throws CorruptSdfException if one of the directories has been deleted since the check
   */
  public static File getRightEnd(final File dir) throws CorruptSdfException {
    if (dir == null || !dir.exists() || !dir.isDirectory()) {
      throw new CorruptSdfException(dir);
    }
    final File rDir = new File(dir, "right");
    if (rDir.exists() && rDir.isDirectory()) {
      return rDir;
    }
    throw new CorruptSdfException(dir);
  }

  /**
   * Checks if the provided reader is empty and <code>NoTalkbackSlimException</code>
   * if it is.
   * @param reader the reader to check.
   * @throws NoTalkbackSlimException if the given sequences reader is null or empty.
   */
  public static void validateNotEmpty(SequencesReader reader) {
    if (reader == null || reader.numberSequences() == 0) {
      final String dir;
      if (reader == null || reader.path() == null) {
        dir = "<Unknown>";
      } else {
        dir = reader.path().getAbsolutePath();
      }
      throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The SDF \"" + dir + "\" was empty.");
    }
  }

  /**
   * Get the GUID from the supplied SDF directory, if the SDF contains one.
   * @param readDir input SDF directory, correctly handles paired end and single end SDF
   * @return GUID for current SDF, 0 if no GUID is specified in the SDF
   * @throws IOException if error reading
   */
  public static SdfId getSdfId(File readDir) throws IOException {
    final File reads = ReaderUtils.isPairedEndDirectory(readDir) ? ReaderUtils.getLeftEnd(readDir) : readDir;
    try (SequencesReader r = SequencesReaderFactory.createDefaultSequencesReader(reads)) {
      return r.getSdfId();
    }
  }

  /**
   * Construct a mapping from sequence name to sequence id.
   * @param sequences the sequences reader to construct the map from
   * @return the map from sequence name to sequence id
   * @throws IOException if an error occurs during reading
   */
  public static Map<String, Long> getSequenceNameMap(final SequencesReader sequences) throws IOException {
    final Map<String, Long> map = new HashMap<>((int) sequences.numberSequences());
    for (long i = 0; i < sequences.numberSequences(); i++) {
      map.put(sequences.name(i), i);
    }
    return map;
  }

  /**
   * Check if a file is an SDF by attempting to get an SDF id
   * @param f file that may be an SDF
   * @return true if the file represents a readable SDF which has an id
   */
  public static boolean isSDF(File f) {
    try {
      return !new SdfId(0).equals(getSdfId(f));
    } catch (IOException e) {
      return false;
    }
  }
}
