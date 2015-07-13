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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(true, false, false), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(false, false, false), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 64, fact.getStateIndex(true, false, false), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 65, fact.getStateIndex(false, false, false), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(true, true, false), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(true, true, true), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(false, true, false), 0, false);
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
    final EvidenceInterface p = fact.evidence(2, 0, 0, 10, 3, fact.getStateIndex(false, true, true), 0, false);
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
    assertEquals(0, fact.getStateIndex(false, false, false));
    assertEquals(1, fact.getStateIndex(false, true, false));
    assertEquals(2, fact.getStateIndex(false, true, true));
    assertEquals(3, fact.getStateIndex(true, false, false));
    assertEquals(4, fact.getStateIndex(true, true, false));
    assertEquals(5, fact.getStateIndex(true, true, true));
  }
}
