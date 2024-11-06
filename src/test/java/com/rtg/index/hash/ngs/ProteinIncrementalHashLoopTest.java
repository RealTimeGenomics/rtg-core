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
package com.rtg.index.hash.ngs;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.io.StringWriter;

import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest.ReadCallMock;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class ProteinIncrementalHashLoopTest extends TestCase {


  public void testReadLoop() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File readDir = new File(dir, "read");
      final File templateDir = new File(dir, "template");
      final String read = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
      final String template = "KWRKNRKSKKNQRNYNHDAADA";
      final int adjLength = read.length() / 3 - 1;
      ReaderTestUtils.getReaderDNA(">a\n" + read, readDir, null).close();
      ReaderTestUtils.getReaderProtein(">b\n" + template, templateDir).close();
      final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new ReadCallMock(new StringWriter()), new ImplementHashFunctionTest.TemplateCallMock());
      final int[] res = new int[3];
      try (ISequenceParams readParams = SequenceParams.builder().directory(readDir).mode(SequenceMode.TRANSLATED).create()) {
        final ProteinIncrementalHashLoop loop = new ProteinIncrementalHashLoop(adjLength, adjLength, mask, false) {
          @Override
          public void hashCall(final int internalId, final int endPosition) {
            res[0]++;
            res[1] += internalId;
            res[2] += endPosition;
          }

          @Override
          public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
            throw new UnsupportedOperationException("Not supported.");
          }


        };
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        final PrintStream ps = new PrintStream(baos);
        Diagnostic.setLogStream(ps);
        try {
          loop.execLoop(readParams);

          ps.flush();
          //          assertTrue(baos.toString().contains("Timer Read_delay"));
        } finally {
          Diagnostic.setLogStream();
        }
      }
      assertEquals(6, res[0]);
      assertEquals(15, res[1]);
      assertEquals(54, res[2]);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testRandomStuff() throws Exception {

    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File readDir = new File(dir, "read");
      final File templateDir = new File(dir, "template");
      final String read = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
      final String template = "KWRKNRKSKKNQRNYNHDAADA";
      final int adjLength = read.length() / 3 - 1;
      ReaderTestUtils.getReaderDNA(">a\n" + read, readDir, null).close();
      ReaderTestUtils.getReaderProtein(">b\n" + template, templateDir).close();
      final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new ReadCallMock(new StringWriter()), new ImplementHashFunctionTest.TemplateCallMock());
      final ProteinIncrementalHashLoop loop = new ProteinIncrementalHashLoop(adjLength, adjLength, mask, false) {

        @Override
        public void hashCall(final int internalId, final int endPosition) {
        }

        @Override
        public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
        }
      };

      try {
        loop.hashCall(3L, 0, 0);
        fail();
      } catch (final UnsupportedOperationException uoe) {
        assertEquals("Don't have a hash to supply here", uoe.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
