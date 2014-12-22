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

import java.io.IOException;
import java.io.InputStream;

import com.rtg.util.io.IOUtils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for <code>Resources</code>.
 *
 */
public class ResourcesTest extends TestCase {

  /**
   * Constructor (needed for JUnit)
   *
   * @param name A string which names the tests.
   *
   */
  public ResourcesTest(final String name) {
    super(name);
  }

  public void testResources() throws IOException {
    InputStream i = Resources.getResourceAsStream("com/rtg/util/ethwinout.txt");
    try {
      IOUtils.readAll(i);
    } finally {
      i.close();
    }
    i = Resources.getResourceAsStream(getClass(), "com/rtg/util/ethwinout.txt");
    try {
      IOUtils.readAll(i);
    } finally {
      i.close();
    }

  }

  public void testSlash() {
    assertEquals("com/rtg/", Resources.trailingSlash("com/rtg", true));
    assertEquals("com/rtg", Resources.trailingSlash("com/rtg", false));
    assertEquals("com/rtg/", Resources.trailingSlash("com/rtg/", true));
    assertEquals("com/rtg", Resources.trailingSlash("com/rtg/", false));
  }
//
//  Unfortunately most of our packages exist in multiple heirachies, making output from listResources non-fixed
//  public void testList(), URISyntaxException {
//    checkList("com/rtg/util");
//    checkList("com/rtg/util/");
//
//  }
//
//  private void checkList(String resource), IOException {
//    String[] res = Resources.listResources(resource);
//    boolean found = false;
//    for (String f : res) {
//      if (f.endsWith("com/rtg/util/ethwinout.txt")) {
//        found = true;
//      }
//    }
//    assertTrue(found);
//  }

  public static Test suite() {
    return new TestSuite(ResourcesTest.class);
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}

