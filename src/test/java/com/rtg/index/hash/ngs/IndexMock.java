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
package com.rtg.index.hash.ngs;


import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;

import com.rtg.index.Finder;
import com.rtg.index.FinderHashValue;
import com.rtg.index.Index;
import com.rtg.index.SparseFrequencyHistogram;
import com.rtg.util.StringUtils;


class IndexMock implements Index {
  private final OutputStreamWriter mOut;
  protected final int mIndex;

  IndexMock(final OutputStream out, final int index) {
    mOut = new OutputStreamWriter(out);
    mIndex = index;
  }
  @Override
  public void add(final long hash, final long id) {
    try {
      mOut.append("add hash=" + hash + " index=" + mIndex + " id=" + id + StringUtils.LS);
      mOut.flush();
    } catch (final IOException e) {
      ReadCallImplementationTest.fail();
    }
  }

  @Override
  public long getInitialHashes() {
    return 0;
  }

  @Override
  public SparseFrequencyHistogram getSparseFrequencyHistogram() {
    return null;
  }

  @Override
  public long bytes() {
    return 0;
  }
  @Override
  public long numberEntries() {
    return 0;
  }
  @Override
  public long numberHashes() {
    return 0;
  }
  @Override
  public void freeze() {
  }
  @Override
  public String infoString() {
    return null;
  }
  @Override
  public String perfString() {
    return null;
  }
  @Override
  public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
    try {
      mOut.append("search hash=").append(String.valueOf(hash)).append(StringUtils.LS);
    } catch (final IOException e) {
      ReadCallImplementationTest.fail();
    }
  }
  @Override
  public void scan(FinderHashValue finder) throws IllegalStateException {
    //do nothing
  }
  @Override
  public int count(final long hash) {
    return 0;
  }

  @Override
  public long first(long hash) throws IllegalStateException {
    return 0;
  }

  @Override
  public long getHash(final long found) {
    throw new UnsupportedOperationException();
  }

  @Override
  public long getValue(long found) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean contains(final long hash) {
    throw new UnsupportedOperationException();
  }
  @Override
  public void dumpValues(final PrintStream out) {
  }

  @Override
  public int maxHashCount() {
    return 10;
  }

}
