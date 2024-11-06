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
