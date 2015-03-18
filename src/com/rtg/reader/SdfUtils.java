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
import java.io.FileNotFoundException;
import java.io.IOException;

import com.rtg.launcher.SequenceParams;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Utility class for SDF stuff.
 */
public final class SdfUtils {

  private SdfUtils() { }

  static final int MAX_NO_DUP_SEQUENCE_COUNT = 100000;

  /**
   * Validate that an SDF directory has names.
   * @param sdf the SDF directory
   */
  public static void validateHasNames(final File sdf) {
    try {
      final boolean hasNames;
      if (ReaderUtils.isPairedEndDirectory(sdf)) {
        hasNames = hasNames(ReaderUtils.getLeftEnd(sdf)) || hasNames(ReaderUtils.getRightEnd(sdf));
      } else {
        hasNames = hasNames(sdf);
      }
      if (!hasNames) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "SDF: " + sdf + " has no name data");
      }
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Unable to find file: " + e.getMessage() + " part of SDF: " + sdf);
    } catch (final IOException e) {
      throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Unable to read SDF: " + sdf + " (" + e.getMessage() + ")");
    }
  }

  private static boolean hasNames(File sdf) throws IOException {
    final IndexFile id = new IndexFile(sdf);
    return id.hasNames();
  }



  /**
   * Validate that a given SDF input has no duplicate sequence names
   * @param noMax set to true to disregard maximum sequence count limit
   * @param params the SequenceParams for the given SDF input
   */
  public static void validateNoDuplicates(final SequenceParams params, boolean noMax) {
    validateNoDuplicates(params.reader(), params.directory(), noMax);
  }

  /**
   * Validate that a given SDF input has no duplicate sequence names
   * @param reader the reader for the SDF
   * @param noMax set to true to disregard maximum sequence count limit
   * @param sdf the SDF file
   */
  public static void validateNoDuplicates(final SequencesReader reader, final File sdf, boolean noMax) {
    if (!reader.hasNames()) {
      throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "SDF: " + sdf + " has no name data");
    }
    if (!noMax && reader.numberSequences() > MAX_NO_DUP_SEQUENCE_COUNT) {
      Diagnostic.warning("Too many sequences to check for duplicate sequence names.");
      return;
    }
    try {
      if (NameDuplicateDetector.checkSequence(reader, null)) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Duplicate sequence names detected in SDF: " + sdf);
      }
    } catch (final IOException e) {
      throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Unable to read SDF: " + sdf + " (" + e.getMessage() + ")");
    }
  }
}
