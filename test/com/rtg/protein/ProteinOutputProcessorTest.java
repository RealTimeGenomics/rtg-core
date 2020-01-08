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
import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.TranslatedFrame;
import com.rtg.ngs.MapFlags;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.util.Environment;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

/**
 */
public class ProteinOutputProcessorTest extends AbstractNanoTest {

  private static final String TEMPLATE_DNA = "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGAAATTATAACCACGACGCAGCAGACGCAG";
  private static final String TEMPLATE_PROTEIN = TestUtils.dnaToProtein(TEMPLATE_DNA);

  /** test template */
  public static final String TEMPLATE_FASTA = ">templateName" + LS
    + TEMPLATE_PROTEIN + LS;

  private static String toFasta(String prefix, String[] seqs) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < seqs.length; ++i) {
      sb.append(">").append(prefix).append(i).append(LS).append(seqs[i]).append(LS);
    }
    return sb.toString();
  }

  private static final String[] READS_PERFECT = {
    "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA",
    "AATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAA",
    "ATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAG",
    "TGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGA",
    DnaUtils.reverseComplement("AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA"),
    DnaUtils.reverseComplement("AATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAA"),
    DnaUtils.reverseComplement("ATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAG"),
    DnaUtils.reverseComplement("TGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGA"),
    "AAAAAAATCAAAGAAATTATAACCACGACGCAGCAG"};
  static final String READS_FASTA_PERFECT = toFasta("testRead", READS_PERFECT);

  static final String POP_EXP_1_ALIGNMENTS_TSV_ID = "pop-exp1-alignments.tsv";
  static final String POP_EXP_UNMAPPED_TSV_ID = "pop-exp-unmapped.tsv";


  protected static void checkUnmappedNoHeader(String expected, File actual) throws IOException {
    assertEquals(expected, getFileFromMatch(actual, "#read-id"));
  }

  public static String getFileFromMatch(File actual, String headerStr) throws IOException {
    assertTrue(actual.exists());
    String actualStr = FileUtils.fileToString(actual);
    actualStr = actualStr.substring(actualStr.indexOf(headerStr));
    return actualStr;
  }

  /**
   * Test method for {@link com.rtg.protein.ProteinOutputProcessor#process(long, java.lang.String, int, int, int, int)}.
   * @throws IOException if error
   * @throws InvalidParamsException
   */
  public final void testProcess() throws IOException, InvalidParamsException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try (final TestDirectory tmp = new TestDirectory("proteinop")) {
      final NgsParams params = createParams(tmp, true);
      try (final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null)) {
        assertTrue(p instanceof ProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+3", 8, 0, 0, 0);
        //p.process(0, "+2", 13, 0, 0, 0);
        p.finish();
      }
      mNano.check(POP_EXP_1_ALIGNMENTS_TSV_ID, getFileFromMatch(new File(tmp, TABULAR_ALIGNMENTS), "#template-name"));
      mNano.check(POP_EXP_UNMAPPED_TSV_ID, getFileFromMatch(new File(tmp, UNMAPPED_FILE), "#read-id"));
      TestUtils.containsAll(mps.toString(),
          "Fast identity percentage : 0%",
          "Overhang identity percentage : 0%",
          "Total alignments : 2",
          "Read cache hits  : 0",
          "Alignments skipped due to offset start   : 0",
          "Alignments skipped due to fast identity  : 0",
          "Alignments done      : 2",
          "Alignments done twice: 0",
          "Alignments failed due to alignment score : 0",
          "Alignments failed due to percent ID      : 0",
          "Alignments failed due to E score         : 0",
          "Alignments failed due to bit score       : 0",
          "Alignments retained : 2"
      );
    }
  }

  private static NgsParams createParams(File tmp, boolean outputprot) throws IOException {
    final NgsFilterParams filter = NgsFilterParams.builder()
    .outputFilter(OutputFilter.PROTEIN_TOPN)
    .useids(true)
    .errorLimit(MapFlags.MAX_SCORE)
    .preFilterMinScore(0)
    .preFilterMinOverlap(0)
    .create();
    return NgsParams.builder()
      .outputParams(NgsOutputParams.builder().outputDir(tmp).outputUnmapped(true).filterParams(filter).outputProteinSequences(outputprot).create())
      .buildFirstParams(getSequenceParamsDna(READS_FASTA_PERFECT))
      .searchParams(getSequenceParamsProtein(TEMPLATE_FASTA))
      .proteinScoringMatrix(new ProteinScoringMatrix())
      .create();
  }

  /**
   * Test method for {@link com.rtg.protein.ProteinOutputProcessor#process(long, java.lang.String, int, int, int, int)}.
   * @throws IOException if error
   * @throws InvalidParamsException
   */
  public final void testProcessNoProteinOutput() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteinop")) {
      final NgsParams params = createParams(tmp, false);
      try (final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null)) {
        assertTrue(p instanceof ProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
        p.finish();
      }
      final String res = FileUtils.fileToString(new File(tmp, "alignments.tsv"));
      assertTrue(res.contains("#Version" + "\t" + Environment.getVersion()));
      assertTrue(res.contains(ProteinOutputProcessor.MAPX_OUTPUT_VERSION_HEADER + "\t" + ProteinOutputProcessor.MAPX_OUTPUT_VERSION));
      assertTrue(res.contains(ProteinOutputProcessor.MAPX_TEMPLATE_SDF_ID_HEADER));
      assertTrue(res.contains(ProteinOutputProcessor.MAPX_READ_SDF_ID_HEADER));
      mNano.check("pop-no-prot-alignments.tsv", getFileFromMatch(new File(tmp, TABULAR_ALIGNMENTS), "#template-name"));
      mNano.check(POP_EXP_UNMAPPED_TSV_ID, getFileFromMatch(new File(tmp, UNMAPPED_FILE), "#read-id"));
    }
  }

  public void testInternalIdToFrameCoversion() {
    assertEquals(6, ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME.length);
    assertEquals(ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[0], 1);
    assertEquals(ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[1], 2);
    assertEquals(ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[2], 3);
    assertEquals(ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[3], -1);
    assertEquals(ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[4], -2);
    assertEquals(ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[5], -3);

    assertEquals(7, ProteinOutputProcessor.FRAMES_MAPPING.length);
    assertEquals(TranslatedFrame.REVERSE3, ProteinOutputProcessor.FRAMES_MAPPING[ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[5] + 3]);
    assertEquals(TranslatedFrame.REVERSE2, ProteinOutputProcessor.FRAMES_MAPPING[ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[4] + 3]);
    assertEquals(TranslatedFrame.REVERSE1, ProteinOutputProcessor.FRAMES_MAPPING[ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[3] + 3]);
    assertEquals(TranslatedFrame.FORWARD3, ProteinOutputProcessor.FRAMES_MAPPING[ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[2] + 3]);
    assertEquals(TranslatedFrame.FORWARD2, ProteinOutputProcessor.FRAMES_MAPPING[ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[1] + 3]);
    assertEquals(TranslatedFrame.FORWARD1, ProteinOutputProcessor.FRAMES_MAPPING[ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[0] + 3]);
    assertEquals(null, ProteinOutputProcessor.FRAMES_MAPPING[3]);
  }

}
