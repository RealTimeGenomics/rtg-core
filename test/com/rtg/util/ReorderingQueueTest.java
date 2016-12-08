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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ReorderingQueueTest extends TestCase {

  static class SimpleRecord {
    final String mRef;
    final int mPos;
    SimpleRecord(String ref, int pos) {
      mRef = ref;
      mPos = pos;
    }
  }

  static class SimpleComparator implements Comparator<SimpleRecord>, Serializable {
    @Override
    public int compare(SimpleRecord o1, SimpleRecord o2) {
      final int r1 = o1.mRef.compareTo(o2.mRef);
      return (r1 != 0) ? r1 : o1.mPos - o2.mPos;
    }
  }

  static class SimpleReorderingQueue extends ReorderingQueue<SimpleRecord> {
    List<SimpleRecord> mOutput = new ArrayList<>();
    SimpleReorderingQueue() {
      super(5, new SimpleComparator());
    }

    @Override
    protected String getReferenceName(SimpleRecord record) {
      return record.mRef;
    }

    @Override
    protected int getPosition(SimpleRecord record) {
      return record.mPos;
    }

    @Override
    protected void flushRecord(SimpleRecord rec) {
      mOutput.add(rec);
    }

    @Override
    protected void reportReorderingFailure(SimpleRecord rec) {
    }
  }

  public void test1() throws IOException {
    final SimpleReorderingQueue q = new SimpleReorderingQueue();
    q.addRecord(new SimpleRecord("a", 4));
    q.addRecord(new SimpleRecord("a", 3));
    q.addRecord(new SimpleRecord("a", 2));
    q.addRecord(new SimpleRecord("a", 2));
    q.addRecord(new SimpleRecord("a", 1));
    q.addRecord(new SimpleRecord("a", 0));
    q.close();
    assertEquals(5, q.mOutput.size());
    for (int i = 0; i < 5; ++i) {
      assertEquals(i, q.mOutput.get(i).mPos);
    }
  }

  public void test2() throws IOException {
    final SimpleReorderingQueue q = new SimpleReorderingQueue();
    q.addRecord(new SimpleRecord("a", 4));
    q.addRecord(new SimpleRecord("a", 3));
    q.addRecord(new SimpleRecord("a", 2));
    q.addRecord(new SimpleRecord("b", 1));
    q.addRecord(new SimpleRecord("b", 0));
    q.close();
    assertEquals(5, q.mOutput.size());
    for (int i = 0; i < 3; ++i) {
      final SimpleRecord r = q.mOutput.get(i);
      assertEquals("a", r.mRef);
      assertEquals(i + 2, r.mPos);
    }
    for (int i = 0; i < 2; ++i) {
      final SimpleRecord r = q.mOutput.get(i + 3);
      assertEquals("b", r.mRef);
      assertEquals(i, r.mPos);
    }
  }
}
