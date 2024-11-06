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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.assembler.graph.Graph;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesIterator;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.diagnostic.Diagnostic;

/**
*/
class ReadPairSource implements AutoCloseable {
  List<SequencesIterator> mIterators;
  List<SequencesReader> mReaders;
  private int mMinInsertSize = -1;
  private int mMaxInsertSize = -1;

  ReadPairSource(SequencesReader... readers) {
    if (readers.length < 1) {
      throw new IllegalArgumentException("No readers supplied");
    }
    mReaders = new ArrayList<>();
    mReaders.addAll(Arrays.asList(readers));
    mIterators = new ArrayList<>(mReaders.size());
    for (SequencesReader current : mReaders) {
      mIterators.add(current.iterator());
    }
    final long readCount = mReaders.get(0).numberSequences();
    for (SequencesReader current : mReaders) {
      if (current.numberSequences() != readCount) {
        throw new IllegalArgumentException("The readers don't have the same number of reads");
      }
    }
  }
  synchronized List<byte[]> nextFragments() throws IOException {
    //System.out.println(Thread.currentThread());
    final List<byte[]> result = new ArrayList<>(mIterators.size());
    for (SequencesIterator reader : mIterators) {
      if (!reader.nextSequence()) {
        return null;
      }
      final byte[] read = new byte[reader.currentLength()];
      reader.readCurrent(read);
      result.add(read);
    }
    return result;
  }

  void reset() {
    for (SequencesIterator reader : mIterators) {
      reader.reset();
    }
  }
  int numberFragments() {
    return mReaders.size();
  }

  long numberReads() {
    return mReaders.get(0).numberSequences();
  }

  @Override
  public void close() {
    for (SequencesReader reader : mReaders) {
      try {
        reader.close();
      } catch (IOException e) {
        Diagnostic.warning(e.getMessage());
      }
    }
  }

  int minInsertSize() {
    return mMinInsertSize;
  }
  int maxInsertSize() {
    return mMaxInsertSize;
  }

  void setMaxInsertSize(int maxInsertSize) {
    mMaxInsertSize = maxInsertSize;
  }

  void setMinInsertSize(int minInsertSize) {
    mMinInsertSize = minInsertSize;
  }

  static ReadPairSource makeSource(File readDir, LongRange region) throws IOException {
    if (ReaderUtils.isPairedEndDirectory(readDir)) {
      final SequencesReader left = SequencesReaderFactory.createDefaultSequencesReader(ReaderUtils.getLeftEnd(readDir), region);
      final SequencesReader right = SequencesReaderFactory.createDefaultSequencesReader(ReaderUtils.getRightEnd(readDir), region);
      return new ReadPairSource(left, right);
    } else {
      final SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(readDir, region);
      return new ReadPairSource(reader);
    }
  }

  /**
   *
   * @param graph the graph to align against
   * @param maxMismatches allowed mismatches
   * @param traverse pre calculated links between contigs
   * @return an aligner of the appropriate type for the reads
   */
  GraphAligner aligner(Graph graph, IntegerOrPercentage maxMismatches, GraphTraversions traverse) {
    return new GraphAligner(graph, maxMismatches, traverse);
  }

}
