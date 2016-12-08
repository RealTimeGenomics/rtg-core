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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.ProgramState;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class AsyncKmerIterableTest extends TestCase {
  public void test() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final String[][] sequences = {
          new String[] {
              "ACGTT"
              , "AGGT"
          }
          , new String[] {
      }
          , new String[] {
              "TTTT"
          , "TTNTT"
          }
      };
      final List<File> files = new ArrayList<>();
      for (int i = 0; i < sequences.length; ++i) {
        final File dir = new File(tmpDir, "" + i);
        files.add(dir);
        final String[] reads = sequences[i];
        final StringBuilder sb = new StringBuilder();
        for (int j = 0; j < reads.length; ++j) {
          sb.append(">")
              .append(j)
              .append(LS)
              .append(reads[j])
              .append(LS);
        }
        ReaderTestUtils.getReaderDNA(sb.toString(), dir, new SdfId());
      }

      final List<ReadPairSource> sources = new ArrayList<>();
      for (File f : files) {
        sources.add(ReadPairSource.makeSource(f, LongRange.NONE));
      }
      try {
        final String[] expected = {
            "ACGT"
            , "CGTT"
            , "AGGT"
            , "TTTT"
        };
        try (AsyncKmerIterable iterable = new AsyncKmerIterable(sources, StringKmer.factory(), 4)) {
          int i = 0;
          for (Kmer k : iterable) {
            assertEquals(expected[i], k.toString());
            ++i;
          }
          assertEquals(expected.length, i);
        }
      } finally {
        for (ReadPairSource source : sources) {
          source.close();
        }
      }

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
  public void testEmpty() throws IOException {
    try (AsyncKmerIterable iterable = new AsyncKmerIterable(Collections.<ReadPairSource>emptyList(), StringKmer.factory(), 3)) {
    assertFalse(iterable.iterator().hasNext());
    }
  }

  public void testThreadAssasination() throws IOException {
    final ReadPairSource source = new DeadlyReadPairSource();
    try {
      try (AsyncKmerIterable iterable = new AsyncKmerIterable(Collections.singletonList(source), StringKmer.factory(), 33)) {
        for (Kmer anIterable : iterable) {
          assertNotNull(anIterable);
        }
        fail();
      }
    } catch (ProgramState.SlimAbortException e) {
      // Check that SimpleThreadPool.terminate has been called.
      // As part of the exception resolution chain the above exception terminate should be called on the thread pool.
      // to allow other threads to resolve themselves. This will result in the program state being reset.
      ProgramState.checkAbort();
    }

  }

  private static class DeadlyReadPairSource extends ReadPairSource {
    byte[][] mBytes = {{1, 2, 3}, {3, 2, 1}};
    int mI;

    DeadlyReadPairSource() throws IOException {
      super(ReaderTestUtils.getReaderDnaMemory(">a" + LS + "G" + LS));
    }

    @Override
    synchronized List<byte[]> nextFragments() throws IOException {
      if (mI < mBytes.length) {
        return Collections.singletonList(mBytes[mI++]);
      } else {
        throw new IOException("Die!!!!");
      }
    }
  }
}
