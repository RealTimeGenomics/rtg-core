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
package com.rtg.metagenomics;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.metagenomics.SpeciesParams.SpeciesParamsBuilder;
import com.rtg.metagenomics.matrix.Matrix;
import com.rtg.metagenomics.matrix.MatrixSimple;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.taxonomy.Taxonomy;
import com.rtg.usage.UsageMetric;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import Jama.EigenvalueDecomposition;
import junit.framework.TestCase;

/**
 */
public class SpeciesTest extends TestCase {

  @Override
  public void setUp() throws Exception {
    super.setUp();
    GlobalFlags.resetAccessedStatus();
    Diagnostic.setLogStream();
  }

  public void test() throws IOException {
    //see speciestest.xls test for details
    final MemoryPrintStream out = new MemoryPrintStream();

    final Frag[] frags = new Frag[2];
    final ArrayList<Integer> l = new ArrayList<>();
    l.add(0);
    l.add(0);
    l.add(1);
    frags[0] = new Frag(l);
    frags[0].setMultiplicity(3);
    l.clear();
    l.add(1);
    l.add(2);
    frags[1] = new Frag(l);

    final SpeciesMap sm = new SpeciesMap();
    sm.id(2);
    sm.id(3);
    sm.id(4);
    sm.id(5);
    sm.id(6);
    final Taxonomy taxonomy = new Taxonomy();
    taxonomy.addNode(1, -1, "root", "none");
    taxonomy.addNode(2, 1, "0", "none");
    taxonomy.addNode(3, 1, "1", "none");
    taxonomy.addNode(4, 1, "2", "none");
    taxonomy.addNode(5, 1, "3", "none");
    taxonomy.addNode(6, 1, "4", "none");
    final MockReaderParams mrp = new MockReaderParams(200, 5, SequenceMode.UNIDIRECTIONAL.codeType()) {
      @Override
      public int[] lengths() {
        return new int[] {4, 5, 6, 7, 0};
      }
    };
    final MemoryPrintStream log = new MemoryPrintStream();
    Diagnostic.setLogStream(log.printStream());
    final SpeciesParams params = new SpeciesParamsBuilder().genome(mrp).minIter(10).create();

    final BlockInfo blockInfo = new BlockInfo(42, null, frags, sm, new long[]{4, 5, 6, 7, 0}, params.verbose());
    final String fragInfo = ""
        + "block 42" + LS
        + "0: 2 {2#3, }" + LS
        + "1: 3 {2#4, }" + LS
        + "2: 4 {2#1, }" + LS
        + "3: 5 {}" + LS
        + "4: 6 {}" + LS
        ;
    assertEquals(fragInfo, blockInfo.fragInfo());
    final Species sp = new Species(blockInfo);
    final BlockResult result = sp.solve(10);
    final UsageMetric usageMetric = new UsageMetric();
    final SpeciesTask t = new SpeciesTask(params, null, usageMetric) {
      {
        mTaxonomy = taxonomy;
        mSpeciesMap = sm;
      }
    };
    //t.result(out.printStream(), result, blockInfo, null, null, false, new SpeciesStatistics(), taxonomy);
    t.result(mrp.reader(), out.lineWriter(), result, blockInfo);
    assertTrue(log.toString(), log.toString().contains("Could not determine length of taxonId \"6\", its either not present or the length is 0" + LS));
    Diagnostic.setLogStream();
    //    System.out.println(out.toString());
    TestUtils.containsAll(out.toString().replace('\t', ' '),
        SpeciesTask.SPECIES_HEADER.replace('\t', ' '),
        "0.6306 0.05727 1.000 0.5771 0.05241 1.000 0.0 0.000 0.000 0 0.00 N 1 2 1 none 0",
        "0.3682 0.01869 1.000 0.4212 0.02138 1.000 0.0 0.000 0.000 0 0.00 N 1 3 1 none 1",
        "0.001189 3.872e-26 1.000 0.001632 5.316e-26 1.000 0.0 0.000 0.000 0 0.00 N 1 4 1 none 2");
    //    System.err.println(log.toString());
  }

  public void test2() throws IOException {
    final File dir = FileUtils.createTempDir("species", "test");
    try {
      Diagnostic.setLogStream();
      //see speciestest.xls test for details
      final MemoryPrintStream out = new MemoryPrintStream();

      final Frag[] frags = new Frag[2];
      final ArrayList<Integer> l = new ArrayList<>();
      l.add(0);
      l.add(0);
      l.add(1);
      frags[0] = new Frag(l);
      frags[0].setMultiplicity(3);
      l.clear();
      l.add(1);
      l.add(1);
      frags[1] = new Frag(l);

      final SpeciesMap sm = new SpeciesMap();
      sm.id(2);
      sm.id(3);
      sm.id(5);

      final Taxonomy taxonomy = new Taxonomy();
      taxonomy.addNode(1, -1, "root", "none");
      taxonomy.addNode(2, 1, "0", "none");
      taxonomy.addNode(3, 1, "1", "none");
      taxonomy.addNode(5, 1, "3", "none");
      final MockReaderParams mrp = new MockReaderParams(200, 4, SequenceMode.UNIDIRECTIONAL.codeType()) {
        @Override
        public int[] lengths() {
          return new int[] {4, 5, 6, 7};
        }
      };
      final SpeciesParams params = new SpeciesParamsBuilder().genome(mrp).minIter(10).create();

      final HashMap<Integer, SpeciesInfo> speciesInfo = new HashMap<>();
      speciesInfo.put(2, new SpeciesInfo(2L, 4L, 4.562));
      speciesInfo.put(3, new SpeciesInfo(7L, 11L, 55));
      final BlockInfo blockInfo = new BlockInfo(42, null, frags, sm, new long[] {4, 11, 7}, params.verbose()); //11 = 5 + 6
      final Species sp = new Species(blockInfo);
      final BlockResult blockResult = sp.solve(10);
      final SpeciesTask t = new SpeciesTask(params, null, null) {
        {
          mTaxonomy = taxonomy;
          mSpeciesInfo.putAll(speciesInfo);
          mSampleFactor = 0.5;
          mSpeciesMap = sm;
        }
      };
      t.result(mrp.reader(), out.lineWriter(), blockResult, blockInfo);
      //    System.out.println(out.toString());
      TestUtils.containsAll(out.toString().replace('\t', ' '),
          SpeciesTask.SPECIES_HEADER.replace('\t',  ' '),
          "0.8621 0.1308 1.000 0.3472 0.05267 1.000 0.0 0.000 0.000 0 4.56 N 1 2 1 none 0",
          "0.1379 0.006867 1.000 0.1528 0.007606 1.000 0.0 0.000 0.000 0 55.00 N 1 3 1 none 1");
      final Matrix hessian = sp.hessianR(blockResult.getR());
      assertEquals(5.3333, hessian.get(0, 0), 1e-4);
      assertEquals(2.6666, hessian.get(0, 1), 1e-4);
      assertEquals(2.6666, hessian.get(1, 0), 1e-4);
      assertEquals(82.3333, hessian.get(1, 1), 1e-4);
      //System.err.println("Hessian:");
      //System.err.println(hessian);
      //Extract eigenvectors and eigenvalues.
      final EigenvalueDecomposition ed = hessian.toJama().eig();
      final Vector r = new Vector(3);
      r.set(0, 1.0);
      r.set(1, 1.0);
      r.set(2, 1.0);
      final Vector stdDev = Species.variance(Species.makeFlatMembership(hessian.rows()), blockInfo, r, ed);
      assertEquals(0.4365, Math.sqrt(stdDev.get(0)), 1e-4);
      assertEquals(0.1111, Math.sqrt(stdDev.get(1)), 1e-4);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testprintAll() throws IOException {
    Diagnostic.setLogStream();
    //see speciestest.xls test for details
    final MemoryPrintStream out = new MemoryPrintStream();

    final Frag[] frags = new Frag[2];
    final ArrayList<Integer> l = new ArrayList<>();
    l.add(0);
    l.add(0);
    l.add(1);
    frags[0] = new Frag(l);
    frags[0].setMultiplicity(3);
    l.clear();
    l.add(1);
    l.add(1);
    frags[1] = new Frag(l);

    final SpeciesMap sm = new SpeciesMap();
    sm.id(2);
    sm.id(3);
    sm.id(4);

    final MockReaderParams mrp = new MockReaderParams(200, 4, SequenceMode.UNIDIRECTIONAL.codeType()) {
      @Override
      public int[] lengths() {
        return new int[] {4, 5, 6, 7};
      }
    };
    final SpeciesParams params = new SpeciesParamsBuilder().genome(mrp).minIter(10).printAll(true).minConfidence(0.0).create();

    final Taxonomy taxonomy = new Taxonomy();
    taxonomy.addNode(1, -1, "root", "none");
    taxonomy.addNode(2, 1, "0", "none");
    taxonomy.addNode(3, 1, "1", "none");
    taxonomy.addNode(4, 1, "3", "none");
    final HashMap<Integer, SpeciesInfo> speciesInfo = new HashMap<>();
    speciesInfo.put(2, new SpeciesInfo(2L, 4L, 4.562));
    speciesInfo.put(3, new SpeciesInfo(7L, 11L, 55));

    final BlockInfo blockInfo = new BlockInfo(42, null, frags, sm, new long[] {4, 11, 7}, params.verbose()); //11 = 5 + 6
    final Species sp = new Species(blockInfo);
    final BlockResult result = sp.solve(10);
    final SpeciesTask t = new SpeciesTask(params, null, null) {
      {
        mTaxonomy = taxonomy;
        mSpeciesInfo.putAll(speciesInfo);
        mSampleFactor = 0.5;
        mSpeciesMap = sm;
      }
    };

    t.result(mrp.reader(), out.lineWriter(), result, blockInfo);
    //    System.out.println(out.toString());
    TestUtils.containsAll(out.toString().replace('\t', ' '),
        "0.8621 0.1308 1.000 0.3472 0.05267 1.000 0.0 0.000 0.000 0 4.56 N 1 2 1 none 0",
        "0.1379 0.006867 1.000 0.1528 0.007606 1.000 0.0 0.000 0.000 0 55.00 N 1 3 1 none 1",
        "0.000 0.000 1.000 0.000 0.000 1.000 0.0 0.000 0.000 0 0.00 N 0 4 1 none 3"
        );
  }

  public void testStdDev() {
    checkStdDev(1.0);
    checkStdDev(1.0, 1.0);
    checkStdDev(1.0, 0.1);
    checkStdDev(1.0, 1.0, 1.0);
    checkStdDev(1.0, 2.0, 0.5);
  }

  private void checkStdDev(final double ... stdDev) {
    //System.err.println("=========================");
    //System.err.println("stdDev=" + IntegralAbstract.toString(stdDev));
    final int n = stdDev.length;
    //construct a Hessian with these std dev.
    final Matrix hessian = new MatrixSimple(n);
    for (int i = 0; i < n; ++i) {
      hessian.set(i, i, 1.0 / Math.pow(stdDev[i], 2));
    }
    checkStdDev(hessian, stdDev);
  }

  private void checkStdDev(final Matrix hessian, final double... stdDev) {
    final int n = hessian.rows();
    //System.err.println("Hessian:");
    //System.err.println(hessian);
    //Extract eigenvectors and eigenvalues.
    final EigenvalueDecomposition ed = hessian.toJama().eig();
    final SpeciesMap sm = new SpeciesMap();
    final long[] genomeLengths = new long[n];
    final Vector r = new Vector(n);
    for (int i = 0; i < n; ++i) {
      sm.id(i);
      genomeLengths[i] = 5;
      r.set(i, 1.0);
    }
    final BlockInfo info = new BlockInfo(42, null, new Frag[0], sm, genomeLengths, false);
    final Vector actual = Species.variance(Species.makeFlatMembership(n), info, r, ed);
    //System.err.println(actual);
    for (int i = 0; i < n; ++i) {
      assertEquals(stdDev[i], Math.sqrt(actual.get(i)), 1e-6);
    }
  }

  public void testStdDevRotated() {
    checkStdDevRotated(1.0, 1.0);
    checkStdDevRotated(1.0, 0.1);
    checkStdDevRotated(0.1, 1.0);
  }

  private void checkStdDevRotated(final double a, final double b) {
    final Matrix hessian = new MatrixSimple(2);
    hessian.set(0, 0, (a + b) / 2.0);
    hessian.set(0, 1, (-a + b) / 2.0);
    hessian.set(1, 0, (-a + b) / 2.0);
    hessian.set(1, 1, (a + b) / 2.0);
    final double x = Math.sqrt(1.0 / (2.0 * a) + 1.0 / (2.0 * b));
    checkStdDev(hessian, x, x);
  }
  static final String BROKEN_TAXONOMY = ""
      + "#RTG taxonomy version 1.0" + LS
      + "#taxID\tparentID\trank\tname" + LS
      + "1\t-1\tno rank\troot" + LS
      + "2\t1\tsuper kingdom\tEukaryota" + LS
      + "3\t4\tsuper kingdom\tBrokenTaxon" + LS
    ;
  static final String BROKEN_TAXONOMY_LOOKUP = ""
      + "3\ta" + LS
    ;

  public void testInvalidTaxonomy() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File species = new File(dir, "species");
      ReaderTestUtils.getReaderDNA(">a" + LS + "ACGT", species, new SdfId());
      final File taxonFile = new File(species, "taxonomy.tsv");
      FileUtils.stringToFile(BROKEN_TAXONOMY, taxonFile);
      try {
        FileUtils.stringToFile(BROKEN_TAXONOMY_LOOKUP, new File(species, "taxonomy_lookup.tsv"));
        final MemoryPrintStream mps = new MemoryPrintStream();
        final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(species);
        final MockSequenceParams mock = new MockSequenceParams(sr, SequenceMode.BIDIRECTIONAL);
        final SpeciesParams params = SpeciesParams.builder().genome(mock.readerParams()).create();
        final SpeciesTask speciesTask = new SpeciesTask(params, mps.outputStream(), new UsageMetric());
        speciesTask.exec();
        fail();
      } catch (NoTalkbackSlimException e) {
        TestUtils.containsAll(e.getMessage(), "The taxonomy in the provided SDF is invalid", "Node 3 does not link to root");

      }
    }
  }
}
