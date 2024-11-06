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

package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.util.StringUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SequencesReaderIteratorTest extends TestCase {
  public void test() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final List<File> files = new ArrayList<>();
      for (int i = 0; i < 4; ++i) {
        final File dir = new File(tmpDir, "" + i);
        files.add(dir);
        ReaderTestUtils.getReaderDNA(">" + i + StringUtils.LS + "AA" + StringUtils.LS, dir, new SdfId());
      }
      final SequencesReaderIterator iterator = new SequencesReaderIterator(files);
      for (int i = 0; i < 4; ++i) {
        assertTrue(iterator.hasNext());
        final SequencesReader reader = iterator.next();
        assertEquals("" + i, reader.name(0));
        reader.close();
      }

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
  public void testEmpty() {
    final SequencesReaderIterator iterator = new SequencesReaderIterator(Collections.emptyList());
    assertFalse(iterator.hasNext());
  }
}
