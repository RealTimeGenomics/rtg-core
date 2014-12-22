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

package com.rtg.reader;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Read names list for right arm that uses left arm names for storage.
 */
public class RightSimplePrereadNames extends SimplePrereadNames {
  private final SimplePrereadNames mLeftNames;

  /**
   * Construct a right arm names holder, using the left arm's names.
   * @param leftNames a {@link SimplePrereadNames} for the left arm
   */
  public RightSimplePrereadNames(SimplePrereadNames leftNames) {
    mLeftNames = leftNames;
  }

  @Override
  public long length() {
    return mLeftNames.length();
  }

  @Override
  public String name(long id) {
    return mLeftNames.name(id);
  }

  @Override
  public void setName(long id, String name) {
    // do nothing - using left read name as storage
  }

  @Override
  public long calcChecksum() {
    return mLeftNames.calcChecksum();
  }

  @Override
  public long bytes() {
    //low balling estimate at 2 pointers and a long per entry. TODO make this more reasonable
    return 0L;
  }

  @Override
  public void writeName(Appendable a, long id) throws IOException {
    a.append(name(id));
  }

  @Override
  public void writeName(OutputStream stream, long id) throws IOException {
    stream.write(mLeftNames.getNameBytes(id));
  }

}
