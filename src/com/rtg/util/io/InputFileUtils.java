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

package com.rtg.util.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Utility class for common checks on lists of input files.
 */
public final class InputFileUtils {

  private InputFileUtils() { }

  /**
   * Method to remove redundant canonical paths from a collection of files.
   * @param files collection of input file names to remove redundant canonical paths from.
   * @return collection of unique files, with first occurrence of a canonical path kept only.
   * @throws IOException if an IOException occurs
   */
  public static List<File> removeRedundantPaths(List<File> files) throws IOException {
    final List<File> out = new ArrayList<>();
    final Set<String> paths = new HashSet<>();
    for (File file : files) {
      if (paths.add(file.getCanonicalPath())) {
        out.add(file);
      }
    }
    return out;
  }

  /**
   * Method to check if two file objects represent the same canonical file.
   * @param f1 first file
   * @param f2 second file
   * @return true if the first and second file are the same canonical file.
   * @throws IOException if an IOException occurs
   */
  public static boolean checkIdenticalPaths(File f1, File f2) throws IOException {
    return f1.getCanonicalPath().equals(f2.getCanonicalPath());
  }
}
