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
package com.rtg.bed;

import java.io.File;
import java.io.IOException;

import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

/**
 */
public final class BedUtils {

  /** BED file suffix */
  public static final String BED_SUFFIX = ".bed";

  private BedUtils() { }

  /**
   * @param f file
   * @return true if file has <code>bed</code> or gzipped <code>bed</code> extension
   */
  public static boolean isBedExtension(File f) {
    return f.getName().endsWith(BED_SUFFIX) || f.getName().endsWith(BED_SUFFIX + FileUtils.GZ_SUFFIX);
  }

  /**
   * Create a tabix index for a BED file
   * @param fileToIndex the BED file
   * @throws IOException if there is a problem
   */
  public static void createBedTabixIndex(File fileToIndex) throws IOException {
    try {
      new TabixIndexer(fileToIndex).saveBedIndex();
    } catch (final IllegalArgumentException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + fileToIndex + ": " + e.getMessage());
      throw e;
    } catch (final UnindexableDataException e) {
      Diagnostic.warning("Cannot produce TABIX index for: " + fileToIndex + ": " + e.getMessage());
    }
  }

  /**
   * Create a ReferenceRegions from the specified BED file
   * @param f a pointer to the bed file
   * @return a new <code>ReferenceRegions</code> or null if the argument is null
   * @throws java.io.IOException when reading the file fails
   */
  public static ReferenceRegions regions(File f) throws IOException {
    if (f != null) {
      try (BedReader reader = BedReader.openBedReader(f, null)) {
        return regions(reader);
      }
    } else {
      return null;
    }
  }

  /**
   * Create a new instance from the specified BED file
   * @param reader the BED reader
   * @return a new <code>ReferenceRegions</code>
   * @throws java.io.IOException when reading the file fails
   */
  public static ReferenceRegions regions(BedReader reader) throws IOException {
    final ReferenceRegions regions = new ReferenceRegions();
    while (reader.hasNext()) {
      regions.add(reader.next());
    }
    return regions;
  }
}
