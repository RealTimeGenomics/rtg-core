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
package com.rtg.sam;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import com.rtg.util.Resources;
import com.rtg.variant.VariantAlignmentRecord;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class DedupifyingIteratorTest extends TestCase {

  private static class VarRecordIterator implements Iterator<VariantAlignmentRecord> {

    private final Iterator<SAMRecord> mIt;
    public VarRecordIterator(Iterator<SAMRecord> it) {
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
    try (InputStream in = Resources.getResourceAsStream("com/rtg/sam/resources/duplicates.sam")) {
      final InputStream exp = Resources.getResourceAsStream("com/rtg/sam/resources/deduplicated.sam");
      try {
        final SAMFileReader inSam = new SAMFileReader(in);
        final SAMFileReader expSam = new SAMFileReader(exp);
        final DedupifyingIterator<VariantAlignmentRecord> inIt = new DedupifyingIterator<>(new VarRecordIterator(inSam.iterator()));
        final Iterator<SAMRecord> expIt = expSam.iterator();
        while (inIt.hasNext() && expIt.hasNext()) {
          final VariantAlignmentRecord inRec = inIt.next();
          final VariantAlignmentRecord expRec = new VariantAlignmentRecord(expIt.next());
          assertTrue(expRec.toString() + " cf." + inRec.toString(), expRec.valueCompareTo(inRec) == 0);
        }
        assertFalse(inIt.hasNext());
        assertFalse(expIt.hasNext());
        assertEquals(6, inIt.numFiltered());
      } finally {
        exp.close();
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

  /**
  public void testSingleEnd() {
    final InputStream in = Resources.getResourceAsStream("com/rtg/sam/resources/duplicatesSE.sam");
    try {
      final InputStream exp = Resources.getResourceAsStream("com/rtg/sam/resources/deduplicatedSE.sam");
      try {
        final SAMFileReader inSam = new SAMFileReader(in);
        final SAMFileReader expSam = new SAMFileReader(exp);
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
