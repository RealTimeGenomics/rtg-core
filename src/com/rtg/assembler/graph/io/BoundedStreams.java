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
import java.io.OutputStream;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.store.StoreDirectory;

/**
 */
@TestClass("com.rtg.assembler.graph.io.BoundedOutputStreamTest")
class BoundedStreams extends BoundedOutputStream {

  private final StoreDirectory mDir;

  private final String mPrefix;

  private final String mSuffix;

  BoundedStreams(final StoreDirectory dir, final int maxLength, final String prefix, final String suffix) {
    super(maxLength);
    mDir = dir;
    mPrefix = prefix;
    mSuffix = suffix;
  }

  @Override
  OutputStream nextProxy(int count) throws IOException {
    return mDir.child(mPrefix + count + mSuffix).outputStream();
  }
}
