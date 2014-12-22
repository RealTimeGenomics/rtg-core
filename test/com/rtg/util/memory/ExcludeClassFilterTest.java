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
package com.rtg.util.memory;


import junit.framework.TestCase;

/**
 *
 */
public class ExcludeClassFilterTest extends TestCase {

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(ExcludeClassFilterTest.class);
  }

  public void testFilter() {
    final ExcludeClassFilter filter = new ExcludeClassFilter(new Class<?>[] {Object.class, Double.class});
    assertFalse(filter.accept(null, new Object()));
    assertTrue(filter.accept(null, "bob"));
    assertFalse(filter.accept(null, 40.3));
    final ExcludeClassFilter filter2 = new ExcludeClassFilter(null);
    assertTrue(filter2.accept(null, new Object()));
    filter2.setExcludedClasses(new Class<?>[] {String.class});
    assertFalse(filter2.accept(null, "bob"));
    assertTrue(filter2.accept(null, 40.3));
  }
}
