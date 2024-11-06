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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import com.rtg.sam.SamFilterParams.SamFilterParamsBuilder;
import com.rtg.tabix.IndexUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.Pair;
import com.rtg.util.Populator;
import com.rtg.util.QuickSort;
import com.rtg.util.QuickSortIntIntProxy;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;
import com.rtg.variant.VariantTestUtils;

import junit.framework.TestCase;

/**
 * Test class
 */
public class CircularBufferMultifileSinglePassReaderWindowTest extends TestCase {

  /**
   * Utility method to create a RecordIterator, for testing purposes.
   *
   * @param <V> the record type
   * @param samFiles files to read records from
   * @param filterParams parameters for filtering
   * @param pop the record populator
   * @return a RecordIterator
   * @throws IOException if an IO error occurs
   */
  public static <V extends ReaderRecord<V> & MateInfo> RecordIterator<V> defaultIterator(Collection<File> samFiles, SamFilterParams filterParams, Populator<V> pop) throws IOException {
    final SingletonPopulatorFactory<V> pf = new SingletonPopulatorFactory<>(pop);
    final SamReadingContext context = new SamReadingContext(samFiles, 4, filterParams, SamUtils.getUberHeader(samFiles), null);
    RecordIterator<V> it = new ThreadedMultifileIterator<>(context, pf);
    if (filterParams.findAndRemoveDuplicates()) {
      it = new DedupifyingRecordIterator<>(it);
    }
    return it;
  }

  private static final int[] EXPECTED_FIRST_OVERLAP_5K_8K = {4902, 4928};
  private static final int[] EXPECTED_FIRST_OVERLAP_56K_8K = {5530, 5591, 5596};
  private static final int[] EXPECTED_START_5K_8K_1 = {5184, 5285, 5750, 5881, 5919, 6027, 6141, 6156, 6266, 6286, 6450, 6580, 6607, 6744, 7056, 7175, 7320, 7453};
  private static final int[] EXPECTED_START_5K_8K_2 = {5273, 5391, 5767, 5887, 6048, 6173, 6198, 6317, 6712, 6854, 7706, 7829, 7875, 7996, 7998};
  private static final int[] EXPECTED_START_5K_8K_3 = {5030, 5484, 5591, 5609, 5612, 5733, 5746, 5824, 5934, 6436, 6477, 6480, 6545, 6593, 6622, 7122, 7256};
  private static final int[] EXPECTED_START_5K_8K_4 = {5480, 5530, 5596, 5650, 5763, 5854, 5856, 5890, 5895, 5964, 5980, 5985, 6033, 6069, 6103, 6126, 6230, 6238, 6268, 6394, 6621, 6723, 6743, 6757, 6820, 6848, 6867, 6948, 7490, 7634, 7646, 7696, 7778, 7805, 7806, 7905};

  private static final int[] EXPECTED_START_5K_8K_COMBINED;
  private static final int[] EXPECTED_GENOME_5K_8K;
  static {
    //combine start positions for the 3 files into 1 array
    EXPECTED_START_5K_8K_COMBINED = new int[EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length + EXPECTED_START_5K_8K_3.length + EXPECTED_START_5K_8K_4.length];
    EXPECTED_GENOME_5K_8K = new int[EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length + EXPECTED_START_5K_8K_3.length + EXPECTED_START_5K_8K_4.length];
    System.arraycopy(EXPECTED_START_5K_8K_1, 0, EXPECTED_START_5K_8K_COMBINED, 0, EXPECTED_START_5K_8K_1.length);
    System.arraycopy(EXPECTED_START_5K_8K_2, 0, EXPECTED_START_5K_8K_COMBINED, EXPECTED_START_5K_8K_1.length, EXPECTED_START_5K_8K_2.length);
    System.arraycopy(EXPECTED_START_5K_8K_3, 0, EXPECTED_START_5K_8K_COMBINED, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length, EXPECTED_START_5K_8K_3.length);
    System.arraycopy(EXPECTED_START_5K_8K_4, 0, EXPECTED_START_5K_8K_COMBINED, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length + EXPECTED_START_5K_8K_3.length, EXPECTED_START_5K_8K_4.length);
    //fill in expected genome (file id)
    Arrays.fill(EXPECTED_GENOME_5K_8K, 0, EXPECTED_START_5K_8K_1.length, 0);
    Arrays.fill(EXPECTED_GENOME_5K_8K, EXPECTED_START_5K_8K_1.length, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length, 1);
    Arrays.fill(EXPECTED_GENOME_5K_8K, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length + EXPECTED_START_5K_8K_3.length, 2);
    Arrays.fill(EXPECTED_GENOME_5K_8K, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length + EXPECTED_START_5K_8K_3.length, EXPECTED_START_5K_8K_1.length + EXPECTED_START_5K_8K_2.length + EXPECTED_START_5K_8K_3.length + EXPECTED_START_5K_8K_4.length, 0);
    QuickSort.sort(new QuickSortIntIntProxy(EXPECTED_START_5K_8K_COMBINED, EXPECTED_GENOME_5K_8K));
  }

  private File mDir;


  @Override
  protected void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileUtils.createTempDir("test", "multisinglepass");
  }

  @Override
  protected void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> getBuffer(final String name, final int end) throws Exception {
    final File samFile = VariantTestUtils.bgzipAndIndexResource("com/rtg/sam/resources/" + name + ".sam", mDir);
    final File[] samFiles = {samFile};
    final List<File> list = Arrays.asList(samFiles);
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0, "a", "b", "c");
    final SamRegionRestriction restriction = new SamRegionRestriction("simulatedSequence1", 0, end);
    final RecordIterator<VariantAlignmentRecord> it = defaultIterator(list, new SamFilterParamsBuilder().restriction(restriction).create(), pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = new CircularBufferMultifileSinglePassReaderWindow<>(it, pop, SamUtils.getUberHeader(list).getSequenceIndex("simulatedSequence1"), restriction.getStart(), Integer.MAX_VALUE);
    return new Pair<>(ssrw, it);
  }

  protected Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> getCircularBuffer(final File[] samFiles, final int start, final int end) throws IOException {
    final List<File> list = Arrays.asList(samFiles);
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0, "a", "b", "c");
    final SamRegionRestriction restriction = new SamRegionRestriction("simulatedSequence2", start, end);
    final RecordIterator<VariantAlignmentRecord> it = defaultIterator(list, new SamFilterParamsBuilder().restriction(restriction).create(), pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = new CircularBufferMultifileSinglePassReaderWindow<>(it, pop, SamUtils.getUberHeader(list).getSequenceIndex("simulatedSequence2"), restriction.getStart(), Integer.MAX_VALUE);
    return new Pair<>(ssrw, it);
  }

  public void testSomeMethod() throws IOException {
    final File[] samFiles = {new File(mDir, "samFile1.sam.gz"), new File(mDir, "samFile2.sam.gz"), new File(mDir, "samFile3.sam.gz"), new File(mDir, "samFile4.sam.gz")};
    final File[] tbiFiles = {new File(mDir, "samFile1.sam.gz.tbi"), new File(mDir, "samFile2.sam.gz.tbi"), new File(mDir, "samFile3.sam.gz.tbi"), new File(mDir, "samFile4.sam.gz.tbi")};
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz", samFiles[0]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz", samFiles[1]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow3.sam.gz", samFiles[2]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindowX4.sam.gz", samFiles[3]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz.tbi", tbiFiles[0]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz.tbi", tbiFiles[1]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow3.sam.gz.tbi", tbiFiles[2]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindowX4.sam.gz.tbi", tbiFiles[3]);
    final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> circularBufferPair = getCircularBuffer(samFiles, 5000, 8000);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = circularBufferPair.getA();
    Iterator<VariantAlignmentRecord> it;
    try {
      it = ssrw.recordsOverlap(5000, 8000); // The overlap form of the same query should pick up a couple more than before
      checkStartCoords(it, false, EXPECTED_FIRST_OVERLAP_5K_8K);
      for (int i = 0; i < EXPECTED_START_5K_8K_COMBINED.length; ++i) {
        assertTrue(it.hasNext());
        final ReaderRecord<?> srr = it.next();
        assertEquals(EXPECTED_START_5K_8K_COMBINED[i], srr.getStart() + 1);
        assertEquals(EXPECTED_GENOME_5K_8K[i], srr.getGenome());
      }
      assertFalse(it.hasNext());

      it = ssrw.recordsOverlap(5600, 5601); // The overlap form of the same query should pick up a couple more than before
      checkStartCoords(it, true, EXPECTED_FIRST_OVERLAP_56K_8K);
      assertFalse(it.hasNext());

    } finally {
      ssrw.flush(5000, 8000);
      ssrw.close();
      circularBufferPair.getB().close();
    }
  }

  private void checkStartCoords(Iterator<VariantAlignmentRecord> it, boolean precise, int... startCoords) { // startCoords here are one-based
    for (final int startPos : startCoords) {
      assertTrue(it.hasNext());
      final ReaderRecord<?> srr = it.next();
      assertEquals(startPos, srr.getStart() + 1);
    }
    if (precise) {
      assertFalse(it.hasNext());
    }

  }

  public void testRecordsOverlap() throws IOException, UnindexableDataException {
    File samFile = FileHelper.resourceToGzFile("com/rtg/sam/resources/readerWindowInsert.sam", new File(mDir, "samFile.sam.gz"));
    samFile = IndexUtils.ensureBlockCompressed(samFile);
    new TabixIndexer(samFile, new File(mDir, "samFile.sam.gz.tbi")).saveSamIndex();
    final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> circularBuffer = getCircularBuffer(new File[]{samFile}, 0, 10000);
    //check that can't ask for same start and end (since it isn't supported)
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> a = circularBuffer.getA();

    //TODO re-add check once we turn on complex read query expansion permanently
//    try {
//      a.recordsOverlap(4, 4);
//      fail();
//    } catch (IllegalArgumentException e) {
//
//    }

    final Iterator<VariantAlignmentRecord> variantAlignmentRecordIterator = a.recordsOverlap(3, 4);
    final int[] expectedStarts = {0, 1, 2, 3};
    int i = 0;
    while (variantAlignmentRecordIterator.hasNext()) {
      final VariantAlignmentRecord next = variantAlignmentRecordIterator.next();
      assertEquals(expectedStarts[i++], next.getStart());
    }
    assertEquals(expectedStarts.length, i);
  }

  public void testCallSequence() throws IOException {
    final File[] samFiles = {new File(mDir, "samFile1.sam.gz"), new File(mDir, "samFile2.sam.gz"), new File(mDir, "samFile3.sam.gz"), new File(mDir, "samFile4.sam.gz")};
    final File[] tbiFiles = {new File(mDir, "samFile1.sam.gz.tbi"), new File(mDir, "samFile2.sam.gz.tbi"), new File(mDir, "samFile3.sam.gz.tbi"), new File(mDir, "samFile4.sam.gz.tbi")};
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz", samFiles[0]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz", samFiles[1]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow3.sam.gz", samFiles[2]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindowX4.sam.gz", samFiles[3]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz.tbi", tbiFiles[0]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz.tbi", tbiFiles[1]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow3.sam.gz.tbi", tbiFiles[2]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindowX4.sam.gz.tbi", tbiFiles[3]);
    final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> circularBufferPair = getCircularBuffer(samFiles, 5000, 8000);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = circularBufferPair.getA();
    Iterator<VariantAlignmentRecord> it;

    final FlushLocus l1 = new FlushLocus(500, 5800);
    final FlushLocus l2 = new FlushLocus(5800, 6000);
    assertTrue(l1.isJoinable(l2));
    assertTrue(l2.isJoinable(l1));

    // To check the coords manually:
    // zgrep simulatedSequence2 readerWindow1.sam.gz readerWindow2.sam.gz readerWindow3.sam.gz readerWindowX4.sam.gz | gawk '$4 >= 5881 && $4 < 6001 {print$0}' | sort -n -k 4

    // Complex calling in first block
    it = ssrw.recordsOverlap(5029, 5060); // Has a record starting at 5029 (0-based), ending outside the range
    checkStartCoords(it, true, 5030);
    it = ssrw.recordsOverlap(5040, 5060); // Has a record starting before and ending after the range
    checkStartCoords(it, true, 5030);

    //ssrw.dumpBuffers();
    // Flush first block
    ssrw.flush(5000, 5060);

    //ssrw.dumpBuffers();
    // Complex calling second block
    it = ssrw.recordsOverlap(5600, 5601); // Has several results
    checkStartCoords(it, true, 5530, 5591, 5596);
    it = ssrw.recordsOverlap(5980, 6000); // Has several results
    checkStartCoords(it, true, /* 5881 doesn't overlap ,*/5887, 5890, 5895, 5919, 5934, 5964, 5980, 5985);

    // Flush second block, test some merging of flush regions
    ssrw.flush(5800, 6000);
    ssrw.flush(5500, 5800);

    // Complex calling in third block
    // Should still be able to access results that overlap from 6000 with start positions < 6000
    it = ssrw.recordsOverlap(6000, 6001); // Has several results
    checkStartCoords(it, true, 5919, 5934, 5964, 5980, 5985);

    // Flush third block
    ssrw.flush(6000, 6001);
    // Complex calling in fourth block
    it = ssrw.recordsOverlap(6500, 6550);
    checkStartCoords(it, true, 6436, 6450, 6477, 6480, 6545);

    // Flush fourth block
    ssrw.flush(6500, 6550);

    // Simple calling in fourth block.. (no such block)
    // Complex calling in third block..

    // Flush third block
    try {
      ssrw.flush(7000, 8000);
      fail("Shouldn't let us flush past what we have read in");
    } catch (final IllegalArgumentException e) {
      // Expected
    }
    ssrw.close();
    circularBufferPair.getB().close();
  }

  public void testLenSaysTheRequirementsChanged() throws Exception {
    final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> bufferPair = getBuffer("readerWindowSmall", 20);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = bufferPair.getA();
    Iterator<VariantAlignmentRecord> it;

    it = ssrw.recordsOverlap(0, 5);
    checkStartCoords(it, true, 1);
    it = ssrw.recordsOverlap(5, 10);
    checkStartCoords(it, true, 7);
    it = ssrw.recordsOverlap(10, 15);
    checkStartCoords(it, true, 11);
    it = ssrw.recordsOverlap(15, 20);
    checkStartCoords(it, true, 17);
    ssrw.flush(15, 20);
    ssrw.flush(10, 15);
    ssrw.flush(5, 10);
    ssrw.flush(0, 5);
    ssrw.close();
    bufferPair.getB().close();
  }

  public void testFlushBeforeRead() throws Exception {
    final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> bufferPair = getBuffer("readerWindowSmall", 20);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = bufferPair.getA();
    Iterator<VariantAlignmentRecord> it;
    assertEquals(0, ssrw.flushedTo());

    it = ssrw.recordsOverlap(5, 10);
    checkStartCoords(it, true, 7);
    ssrw.flush(5, 10);
    assertEquals(0, ssrw.flushedTo());

    it = ssrw.recordsOverlap(0, 5);
    checkStartCoords(it, true, 1);

    it = ssrw.recordsOverlap(10, 15);
    checkStartCoords(it, true, 11);

    it = ssrw.recordsOverlap(15, 20);
    checkStartCoords(it, true, 17);

    ssrw.flush(15, 20);
    ssrw.flush(10, 15);
    ssrw.flush(0, 5);
    assertEquals(20, ssrw.flushedTo());
    ssrw.close();
    bufferPair.getB().close();
  }

  //intermediate chunks have no records - trying to emulate a bug detected in big runs
  public void testFlushBug() throws Exception {
    final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> bufferPair = getBuffer("readerWindowSmallGap", 180);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = bufferPair.getA();
    Iterator<VariantAlignmentRecord> it;
    assertEquals(0, ssrw.flushedTo());

    for (int i = 0; i < 150; i += 15) {
      it = ssrw.recordsOverlap(i, i + 15);
      checkStartCoords(it, true);
    }

    it = ssrw.recordsOverlap(150, 165);
    checkStartCoords(it, true, 151, 157);

    for (int i = 0; i < 165; i += 15) {
      ssrw.flush(i, i + 15);
      assertEquals(i + 15, ssrw.flushedTo());
    }
    ssrw.close();
    bufferPair.getB().close();
  }

  //read and flush in all sorts of crazy ways
  public void testFlushPermute() throws Exception {
    for (int n = 0; n < 20; ++n) {
      //System.err.println("n=" + n);
      final Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> bufferPair = getBuffer("readerWindowSmallGap", 180);
      final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw = bufferPair.getA();
      final RandomArrayList<Integer> ra = new RandomArrayList<>(n);
      for (int i = 0; i < 180; i += 15) {
        ra.add(i + 1);
      }
      while (true) {
        final Integer nx = ra.next();
        if (nx == null) {
          break;
        }
        if (nx > 0) {
          final int start = nx - 1;
          final int end = nx - 1 + 15;
          //System.err.println("recordsOverlap start=" + start + " end=" + end);
          ssrw.recordsOverlap(start, end);
          ra.add(-nx);
        } else {
          final int start = -nx - 1;
          final int end = -nx - 1 + 15;
          //System.err.println("flush start=" + start + " end=" + end);
          ssrw.flush(start, end);
        }
      }
      assertEquals(180, ssrw.flushedTo());
      ssrw.close();
      bufferPair.getB().close();
    }
  }

  public void testAddFlushLocus() {
    checkAddFlushLocus(new FlushLocus[] {new FlushLocus(0, 10)}, new FlushLocus(10, 20), new FlushLocus[] {new FlushLocus(0, 20)});
    checkAddFlushLocus(new FlushLocus[] {new FlushLocus(10, 20)}, new FlushLocus(0, 10), new FlushLocus[] {new FlushLocus(0, 20)});
    checkAddFlushLocus(new FlushLocus[] {new FlushLocus(0, 10), new FlushLocus(20, 30)}, new FlushLocus(10, 20), new FlushLocus[] {new FlushLocus(0, 30)});
    checkAddFlushLocus(new FlushLocus[] {new FlushLocus(0, 10), new FlushLocus(20, 30), new FlushLocus(40, 50)}, new FlushLocus(30, 40), new FlushLocus[] {new FlushLocus(0, 10), new FlushLocus(20, 50)});
    checkAddFlushLocus(new FlushLocus[] {new FlushLocus(0, 10), new FlushLocus(30, 40)}, new FlushLocus(20, 30), new FlushLocus[] {new FlushLocus(0, 10), new FlushLocus(20, 40)});
  }

  private void checkAddFlushLocus(final FlushLocus[] fla, final FlushLocus fl, final FlushLocus[] flb) {
    final List<FlushLocus> flla = makeList(fla);
    final List<FlushLocus> fllb = makeList(flb);
    final List<FlushLocus> afl = CircularBufferMultifileSinglePassReaderWindow.addFlushLocus(flla, fl);
    assertEquals(fllb, afl);
    //System.err.println(afl);
  }

  private List<FlushLocus> makeList(final FlushLocus[] fla) {
    final List<FlushLocus> flla = new LinkedList<>();
    Collections.addAll(flla, fla);
    return flla;
  }

  public void testBug1460Depth0() throws IOException {
    final File[] samFiles = {new File(mDir, "samFile1460.sam.gz")};
    final File[] tbiFiles = {new File(mDir, "samFile1460.sam.gz.tbi")};
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz", samFiles[0]);
    FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow1.sam.gz.tbi", tbiFiles[0]);
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0, "a", "b", "c");
    final SamRegionRestriction region = new SamRegionRestriction("simulatedSequence1", 0, 1000);
    final RecordIterator<VariantAlignmentRecord> it = defaultIterator(Arrays.asList(samFiles), new SamFilterParamsBuilder().restriction(region).create(), pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> buf = new CircularBufferMultifileSinglePassReaderWindow<>(it, pop, 0, region.getStart(), 0);
    try {
      int c = 0;
      for (final Iterator<VariantAlignmentRecord> it2 = buf.recordsOverlap(1, 1000); it2.hasNext();) {
        it2.next();
        ++c;
      }
      assertEquals(8, c);
    } finally {
      buf.close();
      it.close();
    }
  }
}
