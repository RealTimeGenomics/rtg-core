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
