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
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class ReadPairSourceMatePairTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testReads() throws IOException {
    final SequencesReader reader1 = ReaderTestUtils.getReaderDnaMemory(ReadPairSource454Test.lines(">1", "AAACCCTTG"));
    final SequencesReader reader2 = ReaderTestUtils.getReaderDnaMemory(ReadPairSource454Test.lines(">1", "ATACGCTTG"));
    final ReadPairSource source = new ReadPairSourceMatePair(reader1, reader2);
    final List<byte[]> fragments = source.nextFragments();
    assertNotNull(fragments);
    ReadPairSourceTest.listEquals(Arrays.asList(DnaUtils.encodeString("CAAGGGTTT"), DnaUtils.encodeString("CAAGCGTAT")), fragments);
    assertNull(source.nextFragments());
  }
}
