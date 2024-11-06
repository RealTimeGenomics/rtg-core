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
    try (AsyncKmerIterable iterable = new AsyncKmerIterable(Collections.emptyList(), StringKmer.factory(), 3)) {
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
