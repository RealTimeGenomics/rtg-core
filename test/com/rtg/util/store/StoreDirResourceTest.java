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
