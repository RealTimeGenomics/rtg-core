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

package com.rtg.assembler.graph.io;

import java.io.IOException;
import java.util.Collections;

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
    GraphWriter.write(graph0, dir1, "test", Collections.emptySet());
    final Graph graph1 = GraphReader.read(dir1);
    GraphReaderTest.checkGraph(graph1);
  }

  public void testDeletions() throws IOException {
    final MutableGraph write = GraphMapCliTest.makeGraph(1, new String[]{"ACGT", "GGGG", "AATA"}, new long[][]{{1, 2}, {2, 3}, {1, 3}, {1, -3}});
    write.deleteContig(2);
    write.deletePath(3);
    final StoreDirectory store = new StoreDirString();
    GraphWriter.writeWithDeleted(write, store, "testing", Collections.emptySet());
    final Graph read = GraphReader.read(store);
    assertEquals("GGGG", ContigString.contigSequenceString(read.contig(2)));
    assertTrue(read.contigDeleted(2));
    assertFalse(read.contigDeleted(1));
    assertFalse(read.contigDeleted(3));
    assertTrue(read.pathDeleted(3));
    assertFalse(read.pathDeleted(4));
  }
}
