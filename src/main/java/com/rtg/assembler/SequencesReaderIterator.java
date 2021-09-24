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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;

/**
*/
class SequencesReaderIterator implements Iterator<SequencesReader> {
  Iterator<File> mIterator;
  SequencesReaderIterator(List<File> mFiles) {
    mIterator = mFiles.iterator();
  }

  @Override
  public boolean hasNext() {
    return mIterator.hasNext();
  }

  @Override
  public SequencesReader next() {
    try {
      return SequencesReaderFactory.createDefaultSequencesReader(mIterator.next());
    } catch (IOException e) {
      throw new RuntimeException("Wrapped because Iterator doesn't allow throwing", e);
    }
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }
}
