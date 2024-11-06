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
