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
package com.rtg.sam;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import com.rtg.AbstractTest;
import com.rtg.util.Resources;
import com.rtg.variant.VariantAlignmentRecord;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 * Test class
 */
public class DedupifyingIteratorTest extends AbstractTest {

  private static class VarRecordIterator implements Iterator<VariantAlignmentRecord> {

    private final Iterator<SAMRecord> mIt;
    VarRecordIterator(Iterator<SAMRecord> it) {
      mIt = it;
    }

    @Override
    public boolean hasNext() {
      return mIt.hasNext();
    }

    @Override
    public VariantAlignmentRecord next() {
      return new VariantAlignmentRecord(mIt.next());
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
  }

  public void testPairedEnd() throws IOException {
    final boolean skipBug = true; // FYI showing current behaviour, set false to demonstrate bug 1612
    try (InputStream exp = Resources.getResourceAsStream("com/rtg/sam/resources/" + (skipBug ? "deduplicated-current.sam" : "deduplicated.sam"))) {
      try (InputStream in = Resources.getResourceAsStream("com/rtg/sam/resources/duplicates.sam")) {
        final SamReader inSam = SamUtils.makeSamReader(in);
        final SamReader expSam = SamUtils.makeSamReader(exp);
        final DedupifyingIterator<VariantAlignmentRecord> inIt = new DedupifyingIterator<>(new VarRecordIterator(inSam.iterator()));
        final Iterator<SAMRecord> expIt = expSam.iterator();
        while (inIt.hasNext() && expIt.hasNext()) {
          final VariantAlignmentRecord inRec = inIt.next();
          final VariantAlignmentRecord expRec = new VariantAlignmentRecord(expIt.next());
          assertTrue(expRec.toString() + " cf." + inRec.toString(), expRec.valueCompareTo(inRec) == 0);
        }
        if (inIt.hasNext()) {
          fail("Unexpected record: " + inIt.next().toString());
        }
        if (expIt.hasNext()) {
          fail("Filtering stopped before producing expected: " + expIt.next().toString());
        }
        assertEquals(skipBug ? 9 : 6, inIt.numFiltered());
      }
    }
  }

  public void testListSimple() {
    final DedupifyingIterator.DedupifyingList<String> list = new DedupifyingIterator.DedupifyingList<>();
    list.add("a", false);
    list.add("b", true);
    list.add("c", true);
    list.add("d", false);
    list.add("e", false);
    list.add("f", true);
    list.add("g", false);
    list.add("h", true);
    list.add("i", false);
    assertTrue(list.contains("b"));
    assertTrue(list.contains("f"));
    //only checks for ones placed in with true
    assertFalse(list.contains("d"));
    assertFalse(list.contains("g"));
    char c = 'a';
    String i = list.removeFirst();
    while (i != null) {
      assertEquals(Character.toString(c++), i);
      i = list.removeFirst();
    }
    assertTrue(list.isEmpty());
  }

  public void testListDedup() {
    final DedupifyingIterator.DedupifyingList<SimpleWrap> list = new DedupifyingIterator.DedupifyingList<>();
    final SimpleWrap a1 = new SimpleWrap("a");
    final SimpleWrap a2 = new SimpleWrap("a");
    final SimpleWrap a3 = new SimpleWrap("a");
    final SimpleWrap b = new SimpleWrap("b");
    final SimpleWrap c1 = new SimpleWrap("c");
    final SimpleWrap c2 = new SimpleWrap("c");
    final SimpleWrap d = new SimpleWrap("d");
    list.add(a1, true);
    list.add(b, true);
    assertTrue(a1.equals(a2));
    assertTrue(list.contains(a1));
    assertTrue(list.contains(a2));
    list.replace(a2);
    list.add(c1, true);
    list.replace(c2);
    list.replace(a3);
    list.add(d, true);
    final List<SimpleWrap> expOrder = Arrays.asList(b, c2, a3, d);
    final Iterator<SimpleWrap> it = expOrder.iterator();
    SimpleWrap i = list.removeFirst();
    while (i != null) {
      assertTrue(it.hasNext());
      assertTrue(it.next() == i);
      i = list.removeFirst();
    }
    assertFalse(it.hasNext());
  }

  private static final class SimpleWrap {
    private final String mWrapped;
    private SimpleWrap(String wrapped) {
      mWrapped = wrapped;
    }

    @Override
    public boolean equals(Object obj) {
      return obj instanceof SimpleWrap && mWrapped.equals(((SimpleWrap) obj).mWrapped);
    }

    @Override
    public int hashCode() {
      return mWrapped.hashCode();
    }

    @Override
    public String toString() {
      return mWrapped;
    }
  }

  /*
  public void testSingleEnd() {
    final InputStream in = Resources.getResourceAsStream("com/rtg/sam/resources/duplicatesSE.sam");
    try {
      final InputStream exp = Resources.getResourceAsStream("com/rtg/sam/resources/deduplicatedSE.sam");
      try {
        final SamReader inSam = SamUtils.makeSamReader(in);
        final SamReader expSam = SamUtils.makeSamReader(exp);
        final Iterator<VariantAlignmentRecord> inIt = new DedupifyingIterator<>(new VarRecordIterator(inSam.iterator()));
        final Iterator<SAMRecord> expIt = expSam.iterator();
        while (inIt.hasNext() && expIt.hasNext()) {
          final VariantAlignmentRecord inRec = inIt.next();
          final VariantAlignmentRecord expRec = new VariantAlignmentRecord(expIt.next());
          assertTrue(expRec.valueCompareTo(inRec) == 0);
        }
        assertFalse(inIt.hasNext());
        assertFalse(expIt.hasNext());
      } finally {
        exp.close();
      }
    } finally {
      in.close();
    }
  }
  */
}
