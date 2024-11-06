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
