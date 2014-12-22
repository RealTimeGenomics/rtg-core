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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.util.bytecompression.CompressedByteArray;

/**
 * Unit test to try using the <code>FileCompressedInputStream</code> and
 * the <code>FileCompressedOutputStream</code> classes with enough data
 * to break the internal integer casting.
 */
public class FileCompressedStreamRegression extends AbstractFileStreamRegression {

  @Override
  protected long calcLength(int range, long elements) {
    final int bitsPerElement = CompressedByteArray.minBits(range);
    final int bytesPerLong = 8;
    final int elementsPerLong = (8 * bytesPerLong) / bitsPerElement;
    return (elements / elementsPerLong + (elements % elementsPerLong > 0 ? 1 : 0) + 1) * bytesPerLong;
  }

  @Override
  protected OutputStream createOutputStream(File file, int range) throws IOException {
    return new FileCompressedOutputStream(file, range);
  }

  @Override
  protected InputStream createInputStream(File file, int range, long elements, boolean seekable) throws IOException {
    return new FileCompressedInputStream(file, range, elements, seekable);
  }

}
