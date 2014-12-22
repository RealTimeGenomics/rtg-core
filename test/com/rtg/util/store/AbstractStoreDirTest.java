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
import java.util.Iterator;
import java.util.SortedSet;

import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractStoreDirTest extends TestCase {

  /**
   * @param fileDir the directory for file stores, if appropriate
   * @return the store directory interface
   * @throws IOException can happen
   */
  protected abstract StoreDirectory getDirectory(File fileDir) throws IOException;

  public void testDir() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final StoreDirectory dir = getDirectory(tmpDir);
      assertFalse(dir.childExists("child1"));
      assertFalse(dir.childExists("foo"));
      final StoreFile sf1 = dir.child("child1");
      assertEquals("child1", sf1.name());
      sf1.outputStream().close();
      assertTrue(dir.childExists("child1"));
      assertFalse(dir.childExists("foo"));

      final StoreFile sf2 = dir.child("child2");
      assertEquals("child2", sf2.name());
      sf2.outputStream().close();

      final StoreFile sf1a = dir.child("child1");
      assertEquals("child1", sf1a.name());

      final SortedSet<String> children = dir.children();
      assertEquals(2, children.size());
      final Iterator<String> it = children.iterator();
      assertTrue(it.hasNext());
      assertEquals("child1", it.next());
      assertTrue(it.hasNext());
      assertEquals("child2", it.next());
      assertFalse(it.hasNext());
    }
  }
}
