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

import java.io.IOException;
import java.util.SortedSet;

/**
 */
public interface StoreDirectory {

  /**
   * Get a child of the parent directory.
   * @param child name of the child.
   * @return the child as a <code>StoreFile</code>.
   * @throws IOException if accessing child.
   */
  StoreFile child(String child) throws IOException;

  /**
   * Check if the child already exists.
   * @param child to be checked.
   * @return true iff the child exists.
   * @throws IOException when checking.
   */
  boolean childExists(String child) throws IOException;
  /**
   * @return all the children of this parent directory.
   * @throws IOException if problems when accessing children.
   */
  SortedSet<String> children() throws IOException;

}
