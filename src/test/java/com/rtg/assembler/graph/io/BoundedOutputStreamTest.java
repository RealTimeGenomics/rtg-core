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

package com.rtg.assembler.graph.io;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.util.store.StoreDirString;

import junit.framework.TestCase;

/**
 */
public class BoundedOutputStreamTest extends TestCase {

  public void test() throws IOException {
    final StoreDirString dir = new StoreDirString();
    final OutputStream os = new BoundedStreams(dir, 10, "test", "");
    final PrintStream ps = new PrintStream(os);

    ps.print("1234");
    ps.print("5678");
    ps.flush();
    ps.print("ABCD");
    ps.print("EFGH");
    ps.flush();
    ps.print("1234");
    ps.print("5678");
    ps.flush();
    ps.close();
    assertEquals("12345678", dir.child("test1").content());
    assertEquals("ABCDEFGH", dir.child("test2").content());
    assertEquals("12345678", dir.child("test3").content());
  }

  //Check exception if too much buffered when do a flush()
  public void testTooBig() {
    final StoreDirString dir = new StoreDirString();
    final OutputStream os = new BoundedStreams(dir, 10, "test", "");
    final PrintStream ps = new PrintStream(os);

    ps.print("12345678AB");
    ps.flush();
    ps.print("12345678AB");
    ps.flush();
    ps.print("12345678AB");
    ps.flush();
    try {
      ps.print("12345678ABC");
      ps.close();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("11:10", e.getMessage());
    }
  }

  //Test boundaries on two succesive flushes - tests setting of total size
  public void testTricky() throws IOException {
    final StoreDirString dir = new StoreDirString();
    final OutputStream os = new BoundedStreams(dir, 10, "test", "");

    try (PrintStream ps = new PrintStream(os)) {
      ps.print("12345678AB");
      ps.flush();
      ps.print("12345678A");
      ps.flush();
      ps.print("B");
      ps.flush();
      assertEquals("12345678AB", dir.child("test1").content());
      assertEquals("12345678AB", dir.child("test2").content());
      assertEquals("[test1, test2]", dir.children().toString());
    }
  }
}
