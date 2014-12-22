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
    String[] resources;
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
