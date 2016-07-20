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
    final AsyncReadPool pool = new AsyncReadPool("name", Collections.<ReadPairSource>emptyList());
    pool.terminate();
  }
}
