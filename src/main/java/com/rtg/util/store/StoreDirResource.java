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

import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.util.Resources;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class StoreDirResource extends IntegralAbstract implements StoreDirectory {

  private final String mResource;

  private final Map<String, StoreFile> mFiles = new HashMap<>();

  /**
   * @param resource name of the root directory for the resources.
   */
  public StoreDirResource(String resource) {
    mResource = resource;
  }

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
    final StoreFile newCh = new StoreFileResource(mResource, child);
    mFiles.put(child, newCh);
    return newCh;
  }

  @Override
  public boolean childExists(String child) throws IOException {
    final String name = mResource + "/" + child;
    final InputStream fos = Resources.getResourceAsStream(name);
    if (fos != null) {
      fos.close();
    }
    return fos != null;
  }

  @Override
  public SortedSet<String> children() throws IOException {
    final String parent = mResource + "/";
    final int length = parent.length();
    final String[] resources;
    try {
      resources = Resources.listResources(parent);
    } catch (final URISyntaxException e) {
      throw new IOException(e);
    }
    final SortedSet<String> resSet = new TreeSet<>();
    for (final String res : resources) {
      final int index = res.indexOf(parent);
      final String ss = res.substring(index + length);
      //System.err.println("parent=" + parent + " length=" + length + " res=" + res + " index=" + index + " ss=" + ss);
      resSet.add(ss);
    }
    return resSet;
  }
}
