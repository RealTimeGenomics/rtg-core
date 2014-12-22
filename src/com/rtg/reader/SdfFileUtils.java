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

/**
 * Manages files used for a binary sequences directory
 */
public final class SdfFileUtils {

  private SdfFileUtils() {
  }

  /**  main index file name  */
  public static final String INDEX_FILENAME = "mainIndex";
  static final String SEQUENCE_DATA_FILENAME = "seqdata";
  static final String SEQUENCE_QUALITY_DATA_FILENAME = "qualitydata";
  static final String LABEL_DATA_FILENAME = "namedata";
  static final String LABEL_SUFFIX_DATA_FILENAME = "suffixdata";
  static final String SEQUENCE_POINTER_FILENAME = "seqpointer";
  static final String LABEL_POINTER_FILENAME = "namepointer";
  static final String LABEL_SUFFIX_POINTER_FILENAME = "suffixpointer";

  //We currently don't roll index files
  static final String SEQUENCE_INDEX_FILENAME = "sequenceIndex0";
  static final String LABEL_INDEX_FILENAME = "nameIndex0";
  static final String LABEL_SUFFIX_INDEX_FILENAME = "suffixIndex0";

  static File sequenceDataFile(final File dir, final int fileNo) {
    return new File(dir, SEQUENCE_DATA_FILENAME + fileNo);
  }

  static File qualityDataFile(final File dir, final int fileNo) {
    return new File(dir, SEQUENCE_QUALITY_DATA_FILENAME + fileNo);
  }

  static File sequencePointerFile(final File dir, final int fileNo) {
    return new File(dir, SEQUENCE_POINTER_FILENAME + fileNo);
  }

  /**
   * Returns the sequence index file for the given directory, these don't roll.
   * @param dir Directory containing binary sequences
   * @return file for index
   */
  static File sequenceIndexFile(final File dir) {
    return new File(dir, SEQUENCE_INDEX_FILENAME);
  }

  static File labelDataFile(final File dir, final int fileNo) {
    return new File(dir, LABEL_DATA_FILENAME + fileNo);
  }
  static File labelPointerFile(final File dir, final int fileNo) {
    return new File(dir, LABEL_POINTER_FILENAME + fileNo);
  }
  static File labelIndexFile(final File dir) {
    return new File(dir, LABEL_INDEX_FILENAME);
  }

  static File labelSuffixDataFile(final File dir, final int fileNo) {
    return new File(dir, LABEL_SUFFIX_DATA_FILENAME + fileNo);
  }
  static File labelSuffixPointerFile(final File dir, final int fileNo) {
    return new File(dir, LABEL_SUFFIX_POINTER_FILENAME + fileNo);
  }
  static File labelSuffixIndexFile(final File dir) {
    return new File(dir, LABEL_SUFFIX_INDEX_FILENAME);
  }
}
