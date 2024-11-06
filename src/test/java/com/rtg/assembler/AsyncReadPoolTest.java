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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.AbstractTest;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class AsyncReadPoolTest extends AbstractTest {
  public void testPool() throws IOException {
    final MemoryPrintStream logStream = new MemoryPrintStream();
    Diagnostic.setLogStream(logStream.printStream());
    final String[] fasta1 = {"ACGT", "CCCG", "ACCC"};
    final String[] fasta2 = {"CCGT", "AAAG", "AATT"};

    final SequencesReader reader1 = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta(fasta1));
    final SequencesReader reader2 = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta(fasta2));

    final AsyncReadPool pool = new AsyncReadPool("name", Arrays.asList(new ReadPairSource(reader1), new ReadPairSource(reader2)));

    assertEquals(2, pool.sources().size());
    final List<byte[]> fragments = pool.sources().get(0).nextFragments();
    assertEquals(1, fragments.size());
    assertEquals("ACGT", DnaUtils.bytesToSequenceIncCG(fragments.get(0)));

    final List<byte[]> fragments2 = pool.sources().get(1).nextFragments();
    assertEquals(1, fragments2.size());
    assertEquals("CCGT", DnaUtils.bytesToSequenceIncCG(fragments2.get(0)));
    pool.terminate();
    // Can't test that both jobs start because
    TestUtils.containsAll(logStream.toString()
      , "name-0 started reading"
      , "name-0 started reading"
      , "name-0 finished reading"
      , "name-1 finished reading"
    );
  }

  public void testEmpty() throws IOException {
    final AsyncReadPool pool = new AsyncReadPool("name", Collections.emptyList());
    pool.terminate();
  }
}
