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
package com.rtg.position;

import java.io.IOException;
import java.io.PrintStream;
import java.io.StringWriter;

import com.rtg.index.Finder;
import com.rtg.index.FinderHashValue;
import com.rtg.index.Index;
import com.rtg.util.StringUtils;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
final class MockIndex extends IntegralAbstract implements Index {

  private final Appendable mOut = new StringWriter();

  @Override
  public void add(final long hash, final long id) {
    try {
      mOut.append(String.valueOf(id)).append(StringUtils.LS);
    } catch (final IOException e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  public long bytes() {
    throw new UnsupportedOperationException();
  }

  @Override
  public void freeze() {
    throw new UnsupportedOperationException();
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
  public String infoString() {
    throw new UnsupportedOperationException();
  }

  @Override
  public long numberEntries() {
    throw new UnsupportedOperationException();
  }

  @Override
  public long numberHashes() {
    throw new UnsupportedOperationException();
  }

  @Override
  public String perfString() {
    throw new UnsupportedOperationException();
  }

  @Override
  public void search(final long hash, final Finder finder) {
    //do nothing
  }

  @Override
  public void scan(FinderHashValue finder) throws IOException, IllegalStateException {
    //do nothing
  }

  @Override
  public boolean contains(final long hash) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int count(final long hash) {
    throw new UnsupportedOperationException();
  }

  @Override
  public long first(long hash) throws IllegalStateException {
    return 0;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append(mOut);
  }

  @Override
  public boolean globalIntegrity() {
    return true;
  }

  @Override
  public boolean integrity() {
    return true;
  }
  @Override
  public void dumpValues(final PrintStream out) {
  }

  @Override
  public int maxHashCount() {
    throw new UnsupportedOperationException("Not supported yet.");
  }


}
