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
      for (int i = 0; i < 4; i++) {
        final File dir = new File(tmpDir, "" + i);
        files.add(dir);
        ReaderTestUtils.getReaderDNA(">" + i + StringUtils.LS + "AA" + StringUtils.LS, dir, new SdfId());
      }
      final SequencesReaderIterator iterator = new SequencesReaderIterator(files);
      for (int i = 0; i < 4; i++) {
        assertTrue(iterator.hasNext());
        final SequencesReader reader = iterator.next();
        assertTrue(reader.nextSequence());
        assertEquals("" + i, reader.currentName());
        reader.close();
      }

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
  public void testEmpty() {
    final SequencesReaderIterator iterator = new SequencesReaderIterator(Collections.<File>emptyList());
    assertFalse(iterator.hasNext());
  }
}
