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

package com.rtg.assembler.graph.io;

import java.io.IOException;
import java.util.Collections;
import java.util.UUID;

import com.rtg.assembler.GraphMapCliTest;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.util.store.StoreDirResource;
import com.rtg.util.store.StoreDirString;
import com.rtg.util.store.StoreDirectory;

import junit.framework.TestCase;

/**
 * Test that can write and then read back.
 */
public class CommonIOTest extends TestCase {

  //read, write, read then see if is correct at end
  public void test() throws IOException {
    final StoreDirectory dir0 = new StoreDirResource("com/rtg/assembler/planning/resources");
    final Graph graph0 = GraphReader.read(dir0);
    final StoreDirectory dir1 = new StoreDirString();
    GraphWriter.write(graph0, dir1, "test", Collections.<UUID> emptySet());
    final Graph graph1 = GraphReader.read(dir1);
    GraphReaderTest.checkGraph(graph1);
  }

  public void testDeletions() throws IOException {
    final MutableGraph write = GraphMapCliTest.makeGraph(1, new String[]{"ACGT", "GGGG", "AATA"}, new long[][]{{1, 2}, {2, 3}, {1, 3}, {1, -3}});
    write.deleteContig(2);
    write.deletePath(3);
    final StoreDirectory store = new StoreDirString();
    GraphWriter.writeWithDeleted(write, store, "testing", Collections.<UUID>emptySet());
    final Graph read = GraphReader.read(store);
    assertEquals("GGGG", ContigString.contigSequenceString(read.contig(2)));
    assertTrue(read.contigDeleted(2));
    assertFalse(read.contigDeleted(1));
    assertFalse(read.contigDeleted(3));
    assertTrue(read.pathDeleted(3));
    assertFalse(read.pathDeleted(4));
  }
}
