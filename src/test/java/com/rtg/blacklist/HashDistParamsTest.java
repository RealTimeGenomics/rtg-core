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

package com.rtg.blacklist;

import java.io.File;

import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 *
 */
public class HashDistParamsTest extends TestCase {

  public void test() {
    new TestParams(HashDistParams.class, HashDistParamsBuilder.class)
      .excludeParams("close")
      .excludeTypeCheck("blacklistThreshold")
      .check();
  }

  public void testBuilder() {
    final HashDistParamsBuilder b = HashDistParams.builder();
    assertEquals(b, b.makeBlacklist(true));
    assertEquals(b, b.installBlacklist(true));
    assertEquals(b, b.threshold(20));
    assertEquals(b, b.blacklistThreshold(22));
    assertEquals(b, b.numberThreads(2));
    assertEquals(b, b.hashMapSizeFactor(0.42));
    assertEquals(b, b.directory(new File("output")));
    assertEquals(b, b.self());
    
    final HashDistParams params = b.create();
    assertTrue(params.makeBlacklist());
    assertTrue(params.installBlacklist());
    assertEquals(20, params.threshold());
    assertEquals(22, params.blacklistThreshold());
    assertEquals(2, params.numberThreads());
    assertEquals(0.42, params.hashMapSizeFactor());
    assertEquals(new File("output"), params.directory());
    assertEquals(new File("output", "foo"), params.file("foo"));
  }

}