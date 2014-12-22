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
package com.rtg.index.hash.ngs;


import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;

import com.rtg.index.Finder;
import com.rtg.index.FinderHashValue;
import com.rtg.index.Index;
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
  public void scan(FinderHashValue finder) throws IOException, IllegalStateException {
    //do nothing
  }
  @Override
  public int searchCount(final long hash) {
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
  public long search(final long hash) {
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
