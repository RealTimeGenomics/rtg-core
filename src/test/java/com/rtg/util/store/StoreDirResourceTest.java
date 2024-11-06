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

package com.rtg.util.store;

import java.io.File;
import java.io.IOException;
import java.util.SortedSet;

/**
 */
public class StoreDirResourceTest extends AbstractStoreDirTest {

  @Override
  protected StoreDirectory getDirectory(File fileDir) {
    return new StoreDirResource("com/rtg/assembler/graph/io/resources");
  }

  @Override
  public void testDir() throws IOException {
    final File unusedFile = new File("storeDirResourceFileDir");
    final StoreDirectory dir = getDirectory(unusedFile);
    assertFalse(unusedFile.exists());
    assertTrue(dir.childExists("contig1"));
    assertFalse(dir.childExists("foo"));
    final StoreFile sf1 = dir.child("child1");
    assertEquals("child1", sf1.name());
    try {
      sf1.outputStream();
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }

    final StoreFile sf2 = dir.child("child2");
    assertEquals("child2", sf2.name());

    final StoreFile sf1a = dir.child("child1");
    assertEquals("child1", sf1a.name());

    final SortedSet<String> children = dir.children();
    assertEquals(8, children.size());
    //System.err.println(children);
    assertEquals("[contig1, contigLong, header1, header2, headerLong, manifest, paths1, pathsLong]", children.toString());
  }

}
