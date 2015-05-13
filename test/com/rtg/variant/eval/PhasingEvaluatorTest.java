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
package com.rtg.variant.eval;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.Pair;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class PhasingEvaluatorTest extends TestCase {
  static List<DetectedVariant> makeVariantList(List<String> variants, StringBuilder sb) {
    final List<DetectedVariant> callList = new ArrayList<>();
    for (String s : variants) {
      final String vartab = s.replaceAll(" ", "\t");
      callList.add(new DetectedVariant(VcfReader.vcfLineToRecord(vartab), 0, RocSortValueExtractor.NULL_EXTRACTOR, false));
      sb.append(vartab).append("\n");
    }
    return callList;

  }
  private static class MockVariantSet implements VariantSet {
    boolean mSent = false;
    Map<VariantSetType, List<DetectedVariant>> mMap;
    StringBuilder mBaselineVcf = new StringBuilder(VcfHeader.MINIMAL_HEADER + "\tSAMPLE\n");
    StringBuilder mCallsVcf = new StringBuilder(VcfHeader.MINIMAL_HEADER + "\tSAMPLE\n");
    public MockVariantSet(List<String> base, List<String> calls) {
      mMap = new EnumMap<>(VariantSetType.class);
      mMap.put(VariantSetType.BASELINE, makeVariantList(base, mBaselineVcf));
      mMap.put(VariantSetType.CALLS, makeVariantList(calls, mCallsVcf));
    }
    @Override
    public Pair<String, Map<VariantSetType, List<DetectedVariant>>> nextSet() {
      final Map<VariantSetType, List<DetectedVariant>> map;
      if (!mSent) {
         map = mMap;
        mSent = true;
      } else {
        map = null;
      }
      return new Pair<>("10", map);
    }

    @Override
    public VcfHeader baseLineHeader() {
      throw new UnsupportedOperationException();
    }

    @Override
    public VcfHeader calledHeader() {
      throw new UnsupportedOperationException();
    }

    @Override
    public int getNumberOfSkippedBaselineVariants() {
      return 0;
    }

    @Override
    public int getNumberOfSkippedCalledVariants() {
      return 0;
    }
  }
  public void testPhasing() throws IOException, UnindexableDataException {
    final int expectedUnphasable = 0;
    final int expectedMisPhasings = 1;
    final int expectedCorrect = 0;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
          "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 0|1"
    ), Arrays.asList(
          "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 1|0"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);

  }
  public void testDoublePhasing() throws IOException, UnindexableDataException {
    final int expectedUnphasable = 0;
    final int expectedMisPhasings = 2;
    final int expectedCorrect = 0;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 0|1"
        , "10 16 . G T 0.0 PASS . GT 0|1"
    ), Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 1|0"
        , "10 16 . G T 0.0 PASS . GT 0|1"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);
  }

  public void testPhasingUnphased() throws IOException, UnindexableDataException {
    // Test that phase counting resumes correctly after un phased call
    final int expectedUnphasable = 0;
    final int expectedMisPhasings = 2;
    final int expectedCorrect = 1;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 0|1"
        , "10 16 . G T 0.0 PASS . GT 0|1"
        , "10 19 . G T 0.0 PASS . GT 0|1"
        , "10 22 . G T 0.0 PASS . GT 0|1"
        , "10 25 . G T 0.0 PASS . GT 1|0"
    ), Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 1|0"
        , "10 16 . G T 0.0 PASS . GT 1|0"
        , "10 19 . G T 0.0 PASS . GT 0/1"
        , "10 22 . G T 0.0 PASS . GT 1|0"
        , "10 25 . G T 0.0 PASS . GT 1|0"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);

  }

  public void testPhasingUnphased2() throws IOException, UnindexableDataException {
    // Test that phase change after an un phased call isn't counted
    final int expectedUnphasable = 0;
    final int expectedMisPhasings = 1;
    final int expectedCorrect = 2;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 0|1"
        , "10 16 . G T 0.0 PASS . GT 0|1"
        , "10 19 . G T 0.0 PASS . GT 0|1"
        , "10 22 . G T 0.0 PASS . GT 0|1"
        , "10 25 . G T 0.0 PASS . GT 1|0"
    ), Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 1|0"
        , "10 16 . G T 0.0 PASS . GT 1|0"
        , "10 19 . G T 0.0 PASS . GT 0/1"
        , "10 22 . G T 0.0 PASS . GT 0|1"
        , "10 25 . G T 0.0 PASS . GT 1|0"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);

  }
  public void testCluster() throws IOException, UnindexableDataException {
    final int expectedUnphasable = 0;
    final int expectedMisPhasings = 2;
    final int expectedCorrect = 3;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 0|1"
        , "10 15 . G T 0.0 PASS . GT 0|1"
        , "10 17 . G T 0.0 PASS . GT 0|1"
        , "10 19 . G T 0.0 PASS . GT 0|1"
        , "10 25 . G T 0.0 PASS . GT 1|0"
    ), Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 0|1"
        , "10 13 . G T 0.0 PASS . GT 1|0"
        , "10 15 . G T 0.0 PASS . GT 1|0"
        , "10 17 . G T 0.0 PASS . GT 0|1"
        , "10 19 . G T 0.0 PASS . GT 1|0"
        , "10 25 . G T 0.0 PASS . GT 0|1"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);
  }

  public void testUnphaseable() throws IOException, UnindexableDataException {
    final int expectedUnphasable = 1;
    final int expectedMisPhasings = 0;
    final int expectedCorrect = 0;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 0/1"
    ), Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 1|0"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);
  }

  public void testFalsePositive() throws IOException, UnindexableDataException {
    final int expectedUnphasable = 0;
    final int expectedMisPhasings = 1;
    final int expectedCorrect = 0;
    final MockVariantSet variants = new MockVariantSet(Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 16 . G T 0.0 PASS . GT 0|1"
    ), Arrays.asList(
        "10 9 . G T 0.0 PASS . GT 1|0"
        , "10 13 . G T 0.0 PASS . GT 1|0"
        , "10 16 . G T 0.0 PASS . GT 1|0"
    ));
    checkPhasing(expectedCorrect, expectedUnphasable, expectedMisPhasings, variants);
  }

  private void checkPhasing(int expectedCorrect, int expectedUnphasable, int expectedMisPhasings, MockVariantSet variants) throws IOException, UnindexableDataException {
    final MemoryPrintStream tp = new MemoryPrintStream();
    final MemoryPrintStream fp = new MemoryPrintStream();
    final MemoryPrintStream fn = new MemoryPrintStream();
    try (final TestDirectory dir = new TestDirectory()) {
      final File calls = FileHelper.stringToGzFile(variants.mCallsVcf.toString(), new File(dir, "calls.vcf.gz"));
      new TabixIndexer(calls).saveVcfIndex();
      final File baseline = FileHelper.stringToGzFile(variants.mBaselineVcf.toString(), new File(dir, "baseline.vcf.gz"));
      new TabixIndexer(baseline).saveVcfIndex();
      final EvalSynchronizer sync = new EvalSynchronizer(variants, tp.outputStream(), fp.outputStream(), fn.outputStream(), null, baseline, calls, RocSortOrder.DESCENDING);
      final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(VcfEvalTaskTest.REF);
      final SequenceEvaluator eval = new SequenceEvaluator(sync, Collections.singletonMap("10", 0L), reader);
      eval.run();
      assertEquals("correctphasings: " + sync.getCorrectPhasings() + ", misphasings: " + sync.getMisPhasings() + ", unphaseable: " + sync.getUnphasable(), expectedCorrect, sync.getCorrectPhasings());
      assertEquals("misphasings: " + sync.getMisPhasings() + ", unphaseable: " + sync.getUnphasable(), expectedMisPhasings, sync.getMisPhasings());
      assertEquals(expectedUnphasable, sync.getUnphasable());
    }
  }
}
