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
package com.rtg.variant.bayes.snp;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.bayes.EvidenceInterface;

import junit.framework.TestCase;

/**
 */
public class EvidenceQFactoryTest extends TestCase {

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
  }

  public void test() {
    Diagnostic.setLogStream();
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(true, false, true, false), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(0.501187, p.error(), 1e-6);
    assertTrue(p.isForward());
    assertFalse(p.isReadPaired());
    assertFalse(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testNotForward() {
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(false, false, true, false), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(0.501187, p.error(), 1e-6);
    assertFalse(p.isForward());
    assertFalse(p.isReadPaired());
    assertFalse(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testHighPhred() {
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 64, fact.getStateIndex(true, false, true, false), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(5.01187e-7, p.error(), 1e-11);
    assertTrue(p.isForward());
    assertFalse(p.isReadPaired());
    assertFalse(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testHighPhredNotForward() {
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 65, fact.getStateIndex(false, false, true, false), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(5.01187e-7, p.error(), 1e-11);
    assertFalse(p.isForward());
    assertFalse(p.isReadPaired());
    assertFalse(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testForwardUnmated() {
    Diagnostic.setLogStream();
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(true, true, true, false), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(0.501187, p.error(), 1e-6);
    assertTrue(p.isForward());
    assertTrue(p.isReadPaired());
    assertFalse(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testForwardMated() {
    Diagnostic.setLogStream();
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(true, true, true, true), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(0.501187, p.error(), 1e-6);
    assertTrue(p.isForward());
    assertTrue(p.isReadPaired());
    assertTrue(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testNotForwardUnmated() {
    Diagnostic.setLogStream();
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(false, true, true, false), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(0.501187, p.error(), 1e-6);
    assertFalse(p.isForward());
    assertTrue(p.isReadPaired());
    assertFalse(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testNotForwardMated() {
    Diagnostic.setLogStream();
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(false, true, true, true), 0, false);
    assertTrue(p instanceof EvidenceQ);
    assertEquals(0.1, p.mapError());
    assertEquals(2, p.read());
    assertEquals(0.501187, p.error(), 1e-6);
    assertFalse(p.isForward());
    assertTrue(p.isReadPaired());
    assertTrue(p.isMated());
    assertFalse(p.isUnmapped());
  }

  public void testUnmapped() {
    Diagnostic.setLogStream();
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, 0, 0, true);
    assertTrue(p instanceof EvidenceQ);
    assertTrue(p.isUnmapped());
    assertEquals(1.0, p.mapError());
  }

  public void testGetStateIndex() {
    final CachedEvidenceFactory fact = new EvidenceQFactory();
    assertEquals(0, fact.getStateIndex(false, false, true, false));
    assertEquals(1, fact.getStateIndex(false, true, false, false));
    assertEquals(2, fact.getStateIndex(false, true, false, true));
    assertEquals(3, fact.getStateIndex(false, true, true, false));
    assertEquals(4, fact.getStateIndex(false, true, true, true));
    assertEquals(5, fact.getStateIndex(true, false, true, false));
    assertEquals(6, fact.getStateIndex(true, true, false, false));
    assertEquals(7, fact.getStateIndex(true, true, false, true));
    assertEquals(8, fact.getStateIndex(true, true, true, false));
    assertEquals(9, fact.getStateIndex(true, true, true, true));
  }
}
