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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/**
 * Partition files into groups based on file length.
 */
public final class Partition {

  private Partition() { }

  private static class FileLengthComparator implements Comparator<File>, Serializable {
    @Override
    public int compare(final File f1, final File f2) {
      if (f1 == f2) {
        return 0;
      }
      if (f1.length() > f2.length()) {
        return -1;
      } else if (f1.length() < f2.length()) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  private static final FileLengthComparator FILE_COMPARATOR = new FileLengthComparator();

  private static List<List<File>> binThem(final int binCount, final File... files) {
    // Greedy algorithm - sort by file length, start with biggest, put in emptiest bin
    Arrays.sort(files, FILE_COMPARATOR);
    final ArrayList<List<File>> bins = new ArrayList<>();
    for (int k = 0; k < binCount; k++) {
      bins.add(new ArrayList<File>());
    }
    final long[] usage = new long[binCount];

    // This could be made faster if we kept track of most empty bin some other way
    for (final File f : files) {
      int b = 0;
      long u = usage[0];
      for (int k = 1; k < binCount; k++) {
        if (usage[k] < u) { //this results in 0 length files all getting put into the first empty bin, which MAY result in further empty bins remaining empty.
//        if (usage[k] <= u && bins.get(k).size() < bins.get(b).size()) { // <- as per above, this kind of thing would spread the empty files across empty bins.
          u = usage[k];
          b = k;
        }
      }
      usage[b] += f.length();
      bins.get(b).add(f);
    }

    return bins;
  }

  /**
   * Partition the given files into the specified bins trying to approximately
   * balance the size of the bins according to the lengths of the files.
   * Note that some bins may be empty, in the presence of 0 length files.
   *
   * @param binCount number of bins
   * @param files files to bin
   * @return <code>binCount</code> bins collectively containing all the files
   */
  static List<List<File>> partition(final int binCount, final File... files) {
    // Use our own sorted copy to avoid corrupting someone elses array
    final File[] sort = Arrays.copyOf(files, files.length);
    return binThem(binCount, sort);
  }

  /**
   * Partition the given files into the specified bins trying to approximately
   * balance the size of the bins according to the lengths of the files.
   * Note that some bins may be empty, in the presence of 0 length files.
   *
   * @param binCount number of bins
   * @param files files to bin
   * @return <code>binCount</code> bins collectively containing all the files
   */
  public static List<List<File>> partition(final int binCount, final Collection<File> files) {
    return binThem(binCount, files.toArray(new File[files.size()]));
  }
}
