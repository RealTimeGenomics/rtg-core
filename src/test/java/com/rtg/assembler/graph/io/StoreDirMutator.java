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

package com.rtg.assembler.graph.io;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Set;
import java.util.SortedSet;

import com.rtg.util.Pair;
import com.rtg.util.store.StoreDirectory;
import com.rtg.util.store.StoreFile;
import com.rtg.util.store.StoreFileString;

/**
 * Takes a store directory and repalces text in the contents of the various store files.
 * Handy for tweaking files to test error conditions.
 */
public class StoreDirMutator implements StoreDirectory {

  private final StoreDirectory mProxyDirectory;

  private final Set<Pair<String, String>> mRules;


  /**
   * @param proxyDirectory the store directory that is proxy.
   * @param rules for modifying the text that is read.
   */
  StoreDirMutator(StoreDirectory proxyDirectory, Set<Pair<String, String>> rules) {
    super();
    mProxyDirectory = proxyDirectory;
    mRules = rules;
  }

  @Override
  public StoreFile child(String child) throws IOException {
    final StoreFile childFile = mProxyDirectory.child(child);
    String res = childFile.content();
    for (final Pair<String, String> rule : mRules) {
      res = res.replaceAll(rule.getA(), rule.getB());
    }
    final StoreFileString resChild = new StoreFileString(res);
    try (Writer writer = new OutputStreamWriter(resChild.outputStream())) {
      writer.write(res);
    }
    return resChild;
  }

  @Override
  public boolean childExists(String child) throws IOException {
    return mProxyDirectory.childExists(child);
  }

  @Override
  public SortedSet<String> children() throws IOException {
    return mProxyDirectory.children();
  }

}
