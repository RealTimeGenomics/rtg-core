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

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import com.rtg.util.iterators.Transform;

/**
 *
 */
public class AsyncKmerIterable implements KmerIterable {
  final List<ReadPairSource> mSources;
  private final int mKmerSize;
  private final KmerFactory mFactory;
  private AsyncReadPool mPool = null;

  AsyncKmerIterable(List<ReadPairSource> readFiles, KmerFactory factory, int kmerSize) {
    mSources = readFiles;
    mFactory = factory;
    mKmerSize = kmerSize;
  }

  //The chain of transformers needed
  /*
  private static final File2SequencesReader FILE2SEQUENCESREADER = new File2SequencesReader();
  private static class File2SequencesReader extends Transform<File, SequencesReader> {
    @Override
    public SequencesReader trans(File x) {
      try {
        return SequencesReaderFactory.createDefaultSequencesReader(x);
      } catch (final IOException e) {
        throw new RuntimeException(e);
      }
    }
  }
  */

  private static final SequencesReader2Fragments SEQUENCES_READER2FRAGMENTS = new SequencesReader2Fragments();
  private static class SequencesReader2Fragments extends Transform<AsyncReadSource, Iterator<List<byte[]>>> {
    @Override
    public Iterator<List<byte[]>> trans(AsyncReadSource x) {
      return new ReadIterator(x);
    }
  }
  private static class Fragments2Bytes extends Transform<List<byte[]>, Iterator<byte[]>> {
    @Override
    public Iterator<byte[]> trans(List<byte[]> x) {
      return x.iterator();
    }
  }

  private class Bytes2Kmers extends Transform<byte[], Iterator<Kmer>> {
    @Override
    public Iterator<Kmer> trans(byte[] x) {
      return new KmerIterator(x, mFactory, mKmerSize);
    }
  }



  @Override
  public Iterator<Kmer> iterator() {
    for (ReadPairSource source : mSources) {
      source.reset();
    }
    final AsyncReadPool pool = new AsyncReadPool("IterablePool", mSources);
    mPool = pool;
    final Iterator<List<byte[]>> fragments = Transform.flatten(pool.sources().iterator(), SEQUENCES_READER2FRAGMENTS);
    final Iterator<byte[]> bytes = Transform.flatten(fragments, new Fragments2Bytes());
    return Transform.flatten(bytes, new Bytes2Kmers());
  }

  @Override
  public void close() throws IOException {
    mPool.close();
  }
}
