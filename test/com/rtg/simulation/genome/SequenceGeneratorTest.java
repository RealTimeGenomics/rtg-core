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
package com.rtg.simulation.genome;

import java.io.File;

import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Contains tests for related class
 *
 */
public class SequenceGeneratorTest extends TestCase {

  /**
   * Constructor
   */
  public SequenceGeneratorTest(final String name) {
    super(name);
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();

    suite.addTestSuite(SequenceGeneratorTest.class);
    return suite;
  }

  private File mOutDir = null;

  @Override
  public void setUp() throws Exception {
    mOutDir = FileUtils.createTempDir("sequencegenerator", "main");
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    FileHelper.deleteAll(mOutDir);
    mOutDir = null;
  }

  public void testMain() throws Exception {
    final PortableRandom rand = new PortableRandom(1);
    final int[] lengths = {2, 5};
    final int[] freq = {1, 1, 1, 1};
    final RandomDistribution rd = new RandomDistribution(freq, rand);
    final SequenceGenerator sdata = new SequenceGenerator(rand, rd, lengths, mOutDir);
    final long max = sdata.getSizeLimit();
    assertEquals(1000000000, max);
    sdata.createSequences();
    try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(mOutDir)) {
      //System.err.println("" + dsr.maxLength());
      assertEquals(5, dsr.maxLength());
      assertEquals(2, dsr.minLength());
      final String t = "" + dsr.type();
      assertEquals("DNA", t);
      assertEquals("simulatedSequence1", dsr.name(0));
      assertEquals(2, dsr.numberSequences());
    }
  }
}
