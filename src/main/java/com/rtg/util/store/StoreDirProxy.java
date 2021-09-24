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

package com.rtg.util.store;

import java.io.File;
import java.io.IOException;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.io.FileUtils;

/**
 */
public class StoreDirProxy extends IntegralAbstract implements StoreDirectory {

  private final File mDir;

  /**
   * @param dir parent directory
   * @throws IOException if dir is not a directory.
   */
  public StoreDirProxy(final File dir) throws IOException {
    mDir = dir;
    if (!mDir.isDirectory()) {
      throw new IOException("Not a directory " + dir.getPath());
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mDir.isDirectory());
    //System.err.println(path);
    return true;
  }

  @Override
  public StoreFile child(String child) {
    return new StoreFileProxy(new File(mDir, child), child);
  }

  @Override
  public boolean childExists(String child) {
    return new File(mDir, child).isFile();
  }

  @Override
  public SortedSet<String> children() throws IOException {
    final File[] listFiles = FileUtils.listFiles(mDir);
    final SortedSet<String> names = new TreeSet<>();
    for (final File file : listFiles) {
      if (file.isFile()) {
        names.add(file.getName());
      }
    }
    return names;
  }
}
