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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.IORunnable;
import com.rtg.util.Pair;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 *         Date: 12/12/11
 *         Time: 10:57 AM
 */
public class EvalSynchronizerTest extends TestCase {
  private static class MockVariantSet implements VariantSet {
    int mSetId = 0;

    @Override
    public Pair<String, Map<VariantSetType, List<DetectedVariant>>> nextSet() {
      if (mSetId >= 3) {
        return null;
      }
      mSetId++;
      final HashMap<VariantSetType, List<DetectedVariant>> result = new HashMap<>();
      final List<DetectedVariant> empty = Collections.emptyList();
      result.put(VariantSetType.CALLS, empty);
      result.put(VariantSetType.BASELINE, empty);
      return new Pair<String, Map<VariantSetType, List<DetectedVariant>>>("name" + mSetId, result);
    }

    @Override
    public VcfHeader baseLineHeader() {
      return null;
    }

    @Override
    public VcfHeader calledHeader() {
      return null;
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

  private static final String REC1_1 = "name1 3 . A G 0.0 PASS . GT 1/1".replaceAll(" ", "\t");
  private static final String REC1_2 = REC1_1.replaceAll("name1", "name2");
  private static final String REC2_1 = "name1 4 . C G 0.0 PASS . GT 1/1".replaceAll(" ", "\t");
  private static final String REC2_2 = REC2_1.replaceAll("name1", "name2");
  private static final String REC3_1 = "name1 5 . A G 0.0 PASS . GT 1/1".replaceAll(" ", "\t");
  private static final String REC3_2 = REC3_1.replaceAll("name1", "name2");
  private static final String REC4_1 = "name1 6 . C G 0.0 PASS . GT 1/1".replaceAll(" ", "\t");
  private static final String REC4_2 = REC4_1.replaceAll("name1", "name2");
  private static final String REC5_1 = "name1 7 . A G 0.0 PASS . GT 1/1".replaceAll(" ", "\t");
  private static final String REC5_2 = REC5_1.replaceAll("name1", "name2");
  private static final String REC6_1 = "name1 8 . C G 0.0 PASS . GT 1/1".replaceAll(" ", "\t");
  private static final String REC6_2 = REC6_1.replaceAll("name1", "name2");

  private static final String NAME1_RECS = REC1_1 + "\n" + REC2_1 + "\n" + REC3_1 + "\n" + REC4_1 + "\n" + REC5_1 + "\n" + REC6_1 + "\n";
  private static final String NAME2_RECS = REC1_2 + "\n" + REC2_2 + "\n" + REC3_2 + "\n" + REC4_2 + "\n" + REC5_2 + "\n" + REC6_2 + "\n";

  private static final String FAKE_VCF = VcfHeader.MINIMAL_HEADER + "\tSAMPLE\n" + NAME1_RECS + NAME2_RECS;

  public void testEvalSynchronizer() throws IOException, UnindexableDataException {
    final MemoryPrintStream mp = new MemoryPrintStream();
    Diagnostic.setLogStream(mp.printStream());
    try {
      final ByteArrayOutputStream tp = new ByteArrayOutputStream();
      final ByteArrayOutputStream fp = new ByteArrayOutputStream();
      final ByteArrayOutputStream fn = new ByteArrayOutputStream();
      try (final TestDirectory dir = new TestDirectory()) {
        final File fake = FileHelper.stringToGzFile(FAKE_VCF, new File(dir, "fake.vcf.gz"));
        new TabixIndexer(fake).saveVcfIndex();
        final EvalSynchronizer sync = new EvalSynchronizer(new MockVariantSet(), tp, fp, fn, null, fake, fake, RocSortOrder.DESCENDING);
        final Pair<String, Map<VariantSetType, List<DetectedVariant>>> pair = sync.nextSet();
        final Pair<String, Map<VariantSetType, List<DetectedVariant>>> pair2 = sync.nextSet();
        assertEquals("name1", pair.getA());
        assertEquals("name2", pair2.getA());

        final SimpleThreadPool simpleThreadPool = new SimpleThreadPool(3, "pool", true);
        simpleThreadPool.execute(new IORunnable() {
          @Override
          public void run() throws IOException {
            sync.write("name2", Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC1_2), 0, RocSortValueExtractor.NULL_EXTRACTOR, false))
                , Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC3_2), 0, RocSortValueExtractor.NULL_EXTRACTOR, false))
                , Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC5_2), 0, RocSortValueExtractor.NULL_EXTRACTOR, false)), null);
            sync.addVariants(10);
          }
        });
        simpleThreadPool.execute(new IORunnable() {
          @Override
          public void run() throws IOException {
            sync.write("name1", Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC2_1), 0, RocSortValueExtractor.NULL_EXTRACTOR, false))
                , Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC4_1), 0, RocSortValueExtractor.NULL_EXTRACTOR, false))
                , Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC6_1), 0, RocSortValueExtractor.NULL_EXTRACTOR, false)), null);
            sync.addVariants(22);
          }
        });
        simpleThreadPool.terminate();

        assertEquals(REC2_1 + "\n" + REC1_2 + "\n", tp.toString());
        assertEquals(REC4_1 + "\n" + REC3_2 + "\n", fp.toString());
        assertEquals(REC6_1 + "\n" + REC5_2 + "\n", fn.toString());
        assertEquals(32, sync.mTotalVariants);
        assertEquals("name3", sync.nextSet().getA());
        assertEquals(null, sync.nextSet());
      }
      assertTrue(mp.toString().contains("Number of baseline variants: 32"));
      assertTrue(mp.toString().contains("Number of baseline variants: 22") || mp.toString().contains("Number of baseline variants: 10"));
    } finally {
      Diagnostic.setLogStream();
    }

  }
  public void testAbort() throws IOException, UnindexableDataException {
    final ByteArrayOutputStream tp = new ByteArrayOutputStream();
    final OutputStream fp = TestUtils.getNullOutputStream();
    final OutputStream fn = TestUtils.getNullOutputStream();
    try (final TestDirectory dir = new TestDirectory()) {
      final File fake = FileHelper.stringToGzFile(FAKE_VCF, new File(dir, "fake.vcf.gz"));
      new TabixIndexer(fake).saveVcfIndex();
      final EvalSynchronizer sync = new EvalSynchronizer(new MockVariantSet(), tp, fp, fn, null, fake, fake, RocSortOrder.DESCENDING);
      final Pair<String, Map<VariantSetType, List<DetectedVariant>>> pair = sync.nextSet();
      final Pair<String, Map<VariantSetType, List<DetectedVariant>>> pair2 = sync.nextSet();
      assertEquals("name1", pair.getA());
      assertEquals("name2", pair2.getA());

      final SimpleThreadPool simpleThreadPool = new SimpleThreadPool(3, "pool", true);
      simpleThreadPool.execute(new IORunnable() {
        @Override
        public void run() throws IOException {
          sync.write("name2", Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC1_2), 0, RocSortValueExtractor.NULL_EXTRACTOR, false)), null, null, null);
          fail("Should have aborted in thread");
        }
      });
      simpleThreadPool.execute(new IORunnable() {
        @Override
        public void run() throws IOException {
          throw new IOException();
        }
      });
      try {
        simpleThreadPool.terminate();
      } catch (final IOException e) {

      }
    }
  }



  public void testInterrupt() throws IOException, InterruptedException, UnindexableDataException {
    final ByteArrayOutputStream tp = new ByteArrayOutputStream();
    final OutputStream fp = TestUtils.getNullOutputStream();
    final OutputStream fn = TestUtils.getNullOutputStream();
    try (final TestDirectory dir = new TestDirectory()) {
      final File fake = FileHelper.stringToGzFile(FAKE_VCF, new File(dir, "fake.vcf.gz"));
      new TabixIndexer(fake).saveVcfIndex();
      final EvalSynchronizer sync = new EvalSynchronizer(new MockVariantSet(), tp, fp, fn, null, fake, fake, RocSortOrder.DESCENDING);
      final Pair<String, Map<VariantSetType, List<DetectedVariant>>> pair = sync.nextSet();
      final Pair<String, Map<VariantSetType, List<DetectedVariant>>> pair2 = sync.nextSet();
      assertEquals("name1", pair.getA());
      assertEquals("name2", pair2.getA());


      final Exception[] internalException = new Exception[1];
      final Thread t = new Thread() {
        @Override
        public void run() {
          try {
            sync.write("name2", Arrays.asList(new DetectedVariant(VcfReader.vcfLineToRecord(REC1_2), 0, RocSortValueExtractor.NULL_EXTRACTOR, false)), null, null, null);
          } catch (final IllegalStateException e) {
            internalException[0] = e;
          } catch (final IOException e) {
          }
        }
      };
      t.start();
      t.interrupt();
      t.join();
      assertEquals("Interrupted. Unexpectedly", internalException[0].getMessage());
    }
  }
}
