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
package com.rtg.simulation.variants;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SampleSimulatorTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  private static final String REF = ">ref" + StringUtils.LS
          + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
          + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
          + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca";

  public void testSampleSimulator() throws IOException {
    final File dir = FileUtils.createTempDir("genomemut2_", "test");
    try {
      final SequencesReader sr = ReaderTestUtils.getReaderDnaMemory(REF);
      final byte[] buffr = new byte[(int) sr.maxLength()];
      final int lenr = sr.read(0, buffr);
      final String sref = DnaUtils.bytesToSequenceIncCG(buffr, 0, lenr);

      // Generate variants
      final FixedStepPopulationVariantGenerator fixed = new FixedStepPopulationVariantGenerator(sr, 10, new Mutator("X"), new PortableRandom(10), 0.5);
      final List<PopulationVariantGenerator.PopulationVariant> variants = fixed.generatePopulation();
      final File popVcf = new File(dir, "popVcf.vcf.gz");
      PopulationVariantGenerator.writeAsVcf(popVcf, null, variants, sr);
      final String popVarStr = FileHelper.gzFileToString(popVcf);
      //System.out.println("-- Population Variants --");
      //System.out.println(popVarStr);
      final String[] popVars = TestUtils.splitLines(StringUtils.grepMinusV(popVarStr, "^#"));
      assertEquals(12, popVars.length);
      for (String line : popVars) {
        assertEquals(8, line.split("\t").length);
      }

      // Generate sample w.r.t variants
      final SampleSimulator genomemut = new SampleSimulator(sr, new PortableRandom(42), ReferencePloidy.DIPLOID);
      final File vcfOutFile = new File(dir, "sample_foo.vcf.gz");
      genomemut.mutateIndividual(popVcf, vcfOutFile, "foo", Sex.EITHER);
      String sampleVcf = FileHelper.gzFileToString(vcfOutFile);
      //System.out.println("-- Including sample foo --");
      //System.out.println(sampleVcf);
      sampleVcf = StringUtils.grepMinusV(sampleVcf, "^#");
      final String[] sampleVars = TestUtils.splitLines(sampleVcf);
      assertEquals(12, sampleVars.length);
      for (String line : sampleVars) {
        final String[] cols = line.split("\t");
        assertEquals(10, cols.length);
        assertTrue(cols[cols.length - 1].contains("|")); // Generated genotypes are phased
      }

      // Generate SDF corresponding to the sample
      final File outsdf = new File(dir, "outsdf");
      final SampleReplayer vr = new SampleReplayer(sr);
      vr.replaySample(vcfOutFile, outsdf, "foo");

      final SequencesReader srOut = SequencesReaderFactory.createMemorySequencesReader(outsdf, true, LongRange.NONE);
      final byte[] buff = new byte[(int) srOut.maxLength()];
      /*
      System.out.println("-- Chromosomes for sample foo --");
      for (int i = 0; i < srOut.numberSequences(); i++) {
        int len = srOut.read(i, buff);
        System.out.println("seq: " + srOut.name(i));
        System.out.println(DnaUtils.bytesToSequenceIncCG(buff, 0, len));
      }
      */
      assertEquals(2, srOut.numberSequences());
      int len = srOut.read(0, buff);
      final String s1 = DnaUtils.bytesToSequenceIncCG(buff, 0, len);
      len = srOut.read(1, buff);
      final String s2 = DnaUtils.bytesToSequenceIncCG(buff, 0, len);
      assertFalse(sref.equals(s1));
      assertFalse(sref.equals(s2));
      assertFalse(s1.equals(s2));

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
