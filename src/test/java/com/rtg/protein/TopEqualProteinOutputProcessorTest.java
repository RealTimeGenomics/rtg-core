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
package com.rtg.protein;

import static com.rtg.ngs.OutputFilterTest.getSequenceParamsDna;
import static com.rtg.ngs.OutputFilterTest.getSequenceParamsProtein;
import static com.rtg.protein.ProteinOutputProcessor.TABULAR_ALIGNMENTS;
import static com.rtg.protein.ProteinOutputProcessor.UNMAPPED_FILE;
import static com.rtg.protein.ProteinOutputProcessorTest.POP_EXP_1_ALIGNMENTS_TSV_ID;
import static com.rtg.protein.ProteinOutputProcessorTest.POP_EXP_UNMAPPED_TSV_ID;
import static com.rtg.protein.ProteinOutputProcessorTest.getFileFromMatch;

import java.io.File;
import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.io.TestDirectory;

/**
 * Test class
 */
public class TopEqualProteinOutputProcessorTest extends AbstractNanoTest {

  static final String POP_EXP_2_ALIGNMENTS_TSV_ID = "pop-exp2-alignments.tsv";

  public final void testTopEqualProteinOutputProcessor() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteintope")) {
      final NgsParams params = createParams(tmp, 5);
      try (final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null)) {
        assertTrue(p instanceof TopEqualProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
        p.finish();
      }

      mNano.check(POP_EXP_1_ALIGNMENTS_TSV_ID, getFileFromMatch(new File(tmp, TABULAR_ALIGNMENTS), "#template-name"));
      mNano.check(POP_EXP_UNMAPPED_TSV_ID, getFileFromMatch(new File(tmp, UNMAPPED_FILE), "#read-id"));
    }
  }

  public static NgsParams createParams(final File tmp, final int topn) throws IOException, InvalidParamsException {
    return createParams(tmp, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, topn);
  }

  public static NgsParams createParams(final File tmp, String dna, String protein, final int topn) throws IOException, InvalidParamsException {
    final NgsFilterParams filterParams = NgsFilterParams.builder().topN(topn).create();
    return NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).outputUnmapped(true).filterParams(filterParams).create())
      .buildFirstParams(getSequenceParamsDna(dna))
      .searchParams(getSequenceParamsProtein(protein))
      .proteinScoringMatrix(new ProteinScoringMatrix())
      .create();
  }

  public final void testTopEqualProteinOutputProcessor2() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteintope")) {
      final NgsParams params = createParams(tmp, 2);
      try (final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null)) {
        assertTrue(p instanceof TopEqualProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
        p.finish();
      }
      mNano.check(POP_EXP_2_ALIGNMENTS_TSV_ID, getFileFromMatch(new File(tmp, TABULAR_ALIGNMENTS), "#template-name"));
      mNano.check(POP_EXP_UNMAPPED_TSV_ID, getFileFromMatch(new File(tmp, UNMAPPED_FILE), "#read-id"));
    }
  }

  public final void testCountOverflow() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteintope")) {
      final NgsParams params = createParams(tmp, 5);
      try (final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null)) {
        assertTrue(p instanceof TopEqualProteinOutputProcessor);
        for (int i = 0; i < 257; ++i) {
          p.process(0, "+1", 0, 0, 0, 0);
        }
        p.process(0, "+1", 8, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
        p.finish();
      }
      mNano.check(POP_EXP_2_ALIGNMENTS_TSV_ID, getFileFromMatch(new File(tmp, TABULAR_ALIGNMENTS), "#template-name"));
      mNano.check(POP_EXP_UNMAPPED_TSV_ID, getFileFromMatch(new File(tmp, UNMAPPED_FILE), "#read-id"));
    }
  }

  public void testnvalidTopN() throws IOException, InvalidParamsException  {
    try (final TestDirectory tmp = new TestDirectory("proteintope")) {
      for (int topn : new int[] {0, 251 }) {
        try {
          final NgsParams params = createParams(tmp, topn);
          OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
          fail();
        } catch (final IllegalArgumentException ex) {
          assertTrue(ex.getMessage().contains("topN"));
        }
      }
    }
  }
}
