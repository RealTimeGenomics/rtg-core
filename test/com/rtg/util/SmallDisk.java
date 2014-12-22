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
package com.rtg.util;


import java.io.File;

/**
 * Routines to maintain a small file system.  Useful for testing out of
 * disk space conditions in code. It assumes you have created and mounted
 * an appropriate filesystem using <code>
 * internal/scripts/junit/mount-small-disk.sh</code>.
 *
 */
public final class SmallDisk {

  private SmallDisk() { }

  /**
   * Return the root of the small file system.  If the small file system is
   * not available then null is returned.
   *
   * @return root of file system or null
   */
  public static File root() {
    final String user = System.getProperty("user.name");
    final String rootDir = "syscheck".equals(user)
      ? "/home3/syscheck/junit-small-disk"
      : (System.getProperty("user.home") + "/java/cartesian/junit-small-disk");
    final File root = new File(rootDir);
    return root.exists() ? root : null;
  }

  /**
   * Just for testing the root.
   *
   * @param args ignored
   */
  public static void main(final String[] args) {
    System.out.println(root().getPath());
  }
}


