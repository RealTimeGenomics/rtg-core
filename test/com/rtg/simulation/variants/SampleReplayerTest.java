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

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;

import junit.framework.TestCase;

/**
 */
public class SampleReplayerTest extends TestCase {


  private static final String REF = ">simulatedSequence1\n"
      + "GAGACTCGGATCCCCGCTTTTACCGTCTAAGCACTCAAGCTGGAGATTACCATACTTAGGCTCATGTAGCCACCCGCGCTCGTAAATTCTCGACATTCCGCAGTGGCAGCCCTATCGCCA"
      + "GATCAATCGTCGCTGTGGAACTAGACCCGCCCTACTGTTGTCGCTGACTCCTTACTCCCTCTCAACGTATTCATACAGGCCCTTTTCGAATAGCTGGCGTGCACGATTGGCGTGTACATG"
      + "CTGCCTTTGAGGCACAGTTTCGTAGGGTTTCGTAGCAATGGCTTTGTCTGAAGAGAAGGG\n";

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  // Play a deletion into one haplotype, and a snp into the overlapping area of the other haplotype
  public void testOverlappingVariants() throws IOException, UnindexableDataException {
    try (final TestDirectory dir = new TestDirectory("samplereplayer")) {
      final SequencesReader sr = ReaderTestUtils.getReaderDnaMemory(REF);
      final byte[] buffr = new byte[(int) sr.maxLength()];
      final int lenr = sr.read(0, buffr);
      final String sref = DnaUtils.bytesToSequenceIncCG(buffr, 0, lenr);

      final File invcf = new File(dir, "variants.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/util/resources/samplereplayertest.vcf", invcf);
      new TabixIndexer(invcf).saveVcfIndex();
      // Generate SDF corresponding to the sample
      final File outsdf = new File(dir, "outsdf");
      final SampleReplayer vr = new SampleReplayer(sr);
      vr.replaySample(invcf, outsdf, "sm_mom");

      final SequencesReader srOut = SequencesReaderFactory.createMemorySequencesReader(outsdf, true, LongRange.NONE);
      final byte[] buff = new byte[(int) srOut.maxLength()];

      /*
      System.out.println("-- Ref --\n" + sref);
      System.out.println("-- Chromosomes for sample foo --");
      for (int i = 0; i < srOut.numberSequences(); ++i) {
        int len = srOut.read(i, buff);
        System.out.println("seq: " + srOut.name(i));
        System.out.println(DnaUtils.bytesToSequenceIncCG(buff, 0, len));
      }
     */

      assertEquals("GCGTGTACATGCTGC", sref.substring(229, 244));
      int len = srOut.read(0, buff);
      assertEquals("GCGTGAGCTC", DnaUtils.bytesToSequenceIncCG(buff, 0, len).substring(229, 239));
      len = srOut.read(1, buff);
      assertEquals("GCGTGTAAATGCTC", DnaUtils.bytesToSequenceIncCG(buff, 0, len).substring(229, 243));
    }
  }
}
