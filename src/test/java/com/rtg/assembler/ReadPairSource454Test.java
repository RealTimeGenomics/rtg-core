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

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class ReadPairSource454Test extends TestCase {
  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }
  public void testAligner() throws IOException {
    final SequencesReader reader1 = ReaderTestUtils.getReaderDnaMemory(lines(">1", "AAACCCTTG"));
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(2, new String[]{}, new long[][]{});
    assertTrue(new ReadPairSource454(reader1).aligner(graph, new IntegerOrPercentage(2), new GraphTraversions(graph)) instanceof GraphFlowAligner);
  }
  static String lines(String... lines) {
    final StringBuilder sb = new StringBuilder();
    for (String line : lines) {
      sb.append(line);
      sb.append(LS);
    }
    return sb.toString();
  }
  public void testReads() throws IOException {
    final SequencesReader reader1 = ReaderTestUtils.getReaderDnaMemory(lines(">1", "AAACCCTTG"));
    final SequencesReader reader2 = ReaderTestUtils.getReaderDnaMemory(lines(">1", "ATACGCTTG"));
    final ReadPairSource source = new ReadPairSource454(reader1, reader2);
    final List<byte[]> fragments = source.nextFragments();
    assertNotNull(fragments);
    ReadPairSourceTest.listEquals(Arrays.asList(DnaUtils.encodeString("AAACCCTTG")
        , DnaUtils.encodeString("CAAGCGTAT")
    ), fragments);
    assertNull(source.nextFragments());
  }
}
