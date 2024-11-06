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

package com.rtg.variant.bayes.multisample.family;

import junit.framework.TestCase;

/**
 */
public class UniqueIdTest extends TestCase {

  public void test() {
    final UniqueId ui = new UniqueId(3);
    ui.globalIntegrity();
    assertEquals(0, ui.numberIdsSoFar());
    assertEquals(-1, ui.id(0));
    assertEquals(-1, ui.id(1));
    assertEquals(-1, ui.id(2));

    assertEquals(0, ui.addId(1));
    assertEquals(0, ui.addId(1));
    ui.globalIntegrity();
    assertEquals(1, ui.numberIdsSoFar());
    assertEquals(-1, ui.id(0));
    assertEquals(0, ui.id(1));
    assertEquals(-1, ui.id(2));

    assertEquals(1, ui.addId(0));
    assertEquals(1, ui.addId(0));
    ui.globalIntegrity();
    assertEquals(2, ui.numberIdsSoFar());
    assertEquals(1, ui.id(0));
    assertEquals(0, ui.id(1));
    assertEquals(-1, ui.id(2));

    assertEquals(2, ui.addId(2));
    assertEquals(2, ui.addId(2));
    ui.globalIntegrity();
    assertEquals(3, ui.numberIdsSoFar());
    assertEquals(1, ui.id(0));
    assertEquals(0, ui.id(1));
    assertEquals(2, ui.id(2));

  }
}
