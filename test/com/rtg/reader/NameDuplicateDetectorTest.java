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
package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class NameDuplicateDetectorTest extends TestCase {

  class ArrayDummyNameReader extends DummySequencesReader {
    private final String[] mNames;
    private int mIndex = -1;
    public ArrayDummyNameReader(String[] names) {
      mNames = names;
    }
    @Override
    public String currentName() throws IllegalStateException, IOException {
      return mNames[mIndex];
    }
    @Override
    public String currentFullName() throws IllegalStateException, IOException {
      return currentName();
    }
    @Override
    public void seek(long sequenceId) {
      mIndex = (int) sequenceId;
    }
    @Override
    public boolean nextSequence() {
      mIndex++;
      return mIndex < mNames.length && mIndex >= 0;
    }
    @Override
    public long numberSequences() {
      return mNames.length;
    }
    @Override
    public PrereadNamesInterface names() {
      return new PrereadNamesInterface() {
        @Override
        public void writeName(OutputStream os, long id) { }
        @Override
        public void writeName(Appendable a, long id) { }
        @Override
        public String name(long id) {
          return mNames[(int) id];
        }
        @Override
        public long length() {
          return mNames.length;
        }
        @Override
        public long calcChecksum() {
          return 0;
        }
        @Override
        public long bytes() {
          return 0;
        }
      };
    }
    @Override
    public long currentSequenceId() {
      return mIndex;
    }
    @Override
    public String fullName(long index) {
      return name(index);
    }
    @Override
    public void reset() {
      mIndex = -1;
    }
  }

  class NoSortingDuplicateDetector extends NameDuplicateDetector {
    public NoSortingDuplicateDetector(int size) {
      super(size);
    }
    @Override
    void sort() { }
  }

  public void testCheckDuplicates() throws IOException {
    final File tempDir = FileUtils.createTempDir("nameDuplicateDetector", "test");
    try {
      final String[] names1 = {"hello", "sionara", "bye"};
      final String[] names2 = {"00 ", "/O "};
      final String[] names3 = {"bye", "sionara", "hello"};
      final ArrayDummyNameReader reader1 = new ArrayDummyNameReader(names1);
      final ArrayDummyNameReader reader2 = new ArrayDummyNameReader(names2);
      final ArrayDummyNameReader reader3 = new ArrayDummyNameReader(names3);
      NameDuplicateDetector detector = new NameDuplicateDetector(5);
      detector.addPair(names1[0], 0, 0);
      detector.addPair(names1[1], 1, 0);
      detector.addPair(names1[2], 2, 0);
      detector.addPair(names2[0], 0, 1);
      detector.addPair(names2[1], 1, 1);
      final File outputFile = new File(tempDir, "output.txt");
      detector = new NameDuplicateDetector(3);
      detector.addPair(names1[0], 0, 0);
      detector.addPair(names1[1], 1, 0);
      detector.addPair(names1[2], 2, 0);
      assertFalse(detector.checkSequenceDuplicates(new SequencesReader[] {null}, outputFile));
      assertFalse(outputFile.exists());
      detector = new NameDuplicateDetector(2);
      detector.addPair(names2[0], 0, 0);
      detector.addPair(names2[1], 1, 0);
      assertFalse(detector.checkSequenceDuplicates(new SequencesReader[] {reader2}, outputFile));
      assertFalse(outputFile.exists());
      detector = new NameDuplicateDetector(1);
      detector.addPair(names1[0], 0, 0);
      assertFalse(detector.checkSequenceDuplicates(new SequencesReader[] {null}, outputFile));
      assertFalse(outputFile.exists());
      detector = new NameDuplicateDetector(0);
      assertFalse(detector.checkSequenceDuplicates(new SequencesReader[] {null}, outputFile));
      assertFalse(outputFile.exists());
      detector = new NameDuplicateDetector(8);
      detector.addPair(names1[0], 0, 0);
      detector.addPair(names1[1], 1, 0);
      detector.addPair(names1[2], 2, 0);
      detector.addPair(names2[0], 0, 1);
      detector.addPair(names2[1], 1, 1);
      detector.addPair(names3[0], 0, 2);
      detector.addPair(names3[1], 1, 2);
      detector.addPair(names3[2], 2, 2);
      assertTrue(detector.checkSequenceDuplicates(new SequencesReader[] {reader1, reader2, reader3}, outputFile));
      assertTrue(outputFile.exists());
      final String[] expected = {"bye" + StringUtils.LS, "hello" + StringUtils.LS, "sionara" + StringUtils.LS};
      TestUtils.containsAll(FileUtils.fileToString(outputFile), expected);
      assertTrue(outputFile.delete());
      detector = new NameDuplicateDetector(6);
      detector.addPair(names1[0], 0, 0);
      detector.addPair(names1[1], 1, 0);
      detector.addPair(names1[2], 2, 0);
      detector.addPair(names3[0], 0, 1);
      detector.addPair(names3[1], 1, 1);
      detector.addPair(names3[2], 2, 1);
      assertTrue(detector.checkSequenceDuplicates(new SequencesReader[] {reader1, reader3}, outputFile));
      assertTrue(outputFile.exists());
      TestUtils.containsAll(FileUtils.fileToString(outputFile), expected);
      assertTrue(outputFile.delete());
      detector = new NoSortingDuplicateDetector(6);
      detector.addPair(names1[0], 0, 0);
      detector.addPair(names1[1], 1, 0);
      detector.addPair(names1[2], 2, 0);
      detector.addPair(names1[2], 2, 1);
      detector.addPair(names1[1], 1, 1);
      detector.addPair(names1[0], 0, 1);
      try {
        detector.checkSequenceDuplicates(new SequencesReader[] {reader1, reader1}, outputFile);
        fail();
      } catch (final RuntimeException e) {
        assertEquals("List sorting failed", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testNullDetector() throws IOException {
    final NameDuplicateDetector detector = NameDuplicateDetector.getNullDetector();
    detector.addPair("Doesn't Matter", 0, 1);
    assertFalse(detector.checkPrereadDuplicates(null, null));
    assertEquals(0, detector.mHashes.length());
    assertEquals(0, detector.mIndexIds.length());
    assertEquals(0, detector.mCount);
  }

  public void testCheckSequence() throws IOException {
    final File tempDir = FileUtils.createTempDir("nameDuplicateDetector", "test");
    try {
      final String[] names1 = {"hello", "sionara", "bye", "00 ", "/O "};
      final String[] names2 = {"hello", "bye", "bye", "00 ", "/O "};
      final ArrayDummyNameReader reader1 = new ArrayDummyNameReader(names1);
      final ArrayDummyNameReader reader2 = new ArrayDummyNameReader(names2);
      final File outputFile = new File(tempDir, "output.txt");
      assertFalse(NameDuplicateDetector.checkSequence(reader1, outputFile));
      assertFalse(outputFile.exists());
      assertTrue(NameDuplicateDetector.checkSequence(reader2, outputFile));
      assertTrue(outputFile.exists());
      assertEquals("bye" + StringUtils.LS, FileUtils.fileToString(outputFile));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

}
