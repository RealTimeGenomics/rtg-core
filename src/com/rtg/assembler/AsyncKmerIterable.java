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
