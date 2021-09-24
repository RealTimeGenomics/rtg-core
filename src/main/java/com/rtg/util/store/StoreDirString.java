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

import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class StoreDirString extends IntegralAbstract implements StoreDirectory {

  private final Map<String, StoreFile> mFiles = new HashMap<>();

  @Override
  public boolean integrity() {
    Exam.assertTrue(mFiles != null);
    return true;
  }

  @Override
  public StoreFile child(String child) {
    final StoreFile ch = mFiles.get(child);
    if (ch != null) {
      return ch;
    }
    final StoreFile newCh = new StoreFileString(child);
    mFiles.put(child, newCh);
    return newCh;
  }

  @Override
  public boolean childExists(String child) {
    return mFiles.containsKey(child);
  }

  @Override
  public SortedSet<String> children() {
    return new TreeSet<>(mFiles.keySet());
  }


}
