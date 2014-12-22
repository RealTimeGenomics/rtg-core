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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.TranslatedFrame;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Environment;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test Class
 */
public class ProteinOutputProcessorTest extends TestCase {

  private static final String TEMPLATE_DNA = "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGAAATTATAACCACGACGCAGCAGACGCAG";
  private static final String TEMPLATE_PROTEIN = TestUtils.dnaToProtein(TEMPLATE_DNA);

  /** test template */
  public static final String TEMPLATE_FASTA = ">templateName" + StringUtils.LS
    + TEMPLATE_PROTEIN + StringUtils.LS;

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

  /** test reads */
  public static final String READS_FASTA_PERFECT;
  static final String TB = "\t";
  static final String EXPECTED_1 = "#template-name" + TB + "frame" + TB + "read-id" + TB
    + "template-start" + TB + "template-end" + TB + "template-length" + TB
    + "read-start" + TB + "read-end" + TB + "read-length" + TB
    + "template-protein" + TB + "read-protein" + TB + "alignment" + TB
    + "identical" + TB + "%identical" + TB + "positive" + TB + "%positive" + TB + "mismatches" + TB
    + "raw-score" + TB + "bit-score" + TB + "e-score" + StringUtils.LS

    + "templateName" + TB + "+1" + TB + "0" + TB
    + "1" + TB + "12" + TB + "22" + TB + "1" + TB + "36" + TB + "36" + TB
    + "kwrknrkskknq" + TB + "kwrknrkskknq" + TB + "kwrknrkskknq" + TB
    // matches   +matches   mismatches
    + "12" + TB + "100" + TB + "12" + TB + "100" + TB + "0" + TB
    // alignment-    bit-   e-scores
    + "-67" + TB + "30.4" + TB + "1.5e-8" + StringUtils.LS

    + "templateName" + TB + "+3" + TB + "1" + TB
    + "2" + TB + "12" + TB + "22" + TB + "3" + TB + "35" + TB + "36" + TB
    + "wrknrkskknq" + TB + "wrknrkskknq" + TB + "wrknrkskknq" + TB
    // matches   +matches   mismatches
    + "11" + TB + "100" + TB + "11" + TB + "100" + TB + "0" + TB
    + "-62" + TB + "28.5" + TB + "5.8e-8" +  StringUtils.LS;
  static {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < READS_PERFECT.length; i++) {
      sb.append(">testRead").append(i).append(StringUtils.LS).append(READS_PERFECT[i]).append(StringUtils.LS);
    }
    READS_FASTA_PERFECT = sb.toString();
  }

  static final String EXPECTED_UNMAPPED = "#read-id\treason-unmapped" + LS
    + "2" + LS
    + "3" + LS
    + "4" + LS
    + "5" + LS
    + "6" + LS
    + "7" + LS
    + "8" + LS;

  /**
   * Test method for {@link com.rtg.protein.ProteinOutputProcessor#process(long, java.lang.String, int, int, int, int)}.
   * @throws IOException if error
   * @throws InvalidParamsException
   */
  public final void testProcess() throws IOException, InvalidParamsException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final File tmp = FileUtils.createTempDir("proteinop", "filter");
    try {
      final File input = new File(tmp, "1");
      ReaderTestUtils.getReaderDNA(READS_FASTA_PERFECT, input, null).close();
      final File input2 = new File(tmp, "2");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, input2).close();

      final NgsFilterParams filter = NgsFilterParams.builder()
      .outputFilter(OutputFilter.PROTEIN_TOPN)
      .useids(true)
      .errorLimit(CommonFlags.MAX_SCORE)
      .preFilterMinScore(0)
      .preFilterMinOverlap(0)
      .create();
      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).outputUnmapped(true).filterParams(filter).create())
        .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).loadNames(true).create())
        .searchParams(SequenceParams.builder().directory(input2).useMemReader(true).loadNames(true).create())
        .proteinScoringMatrix(new ProteinScoringMatrix())
        .create();
      final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null);
      try {
        assertTrue(p instanceof ProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
      } finally {
        p.finish();
        p.close();
      }
      MapXCliTest.checkAlignmentsNoHeader(EXPECTED_1, new File(tmp, "alignments.tsv"));
      MapXCliTest.checkUnmappedNoHeader(EXPECTED_UNMAPPED, new File(tmp, ProteinOutputProcessor.UNMAPPED_FILE));
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
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  static final String EXPECTED_NO_PROTEIN = "#template-name" + TB + "frame" + TB + "read-id" + TB
  + "template-start" + TB + "template-end" + TB + "template-length" + TB
  + "read-start" + TB + "read-end" + TB + "read-length" + TB
  + "identical" + TB + "%identical" + TB + "positive" + TB + "%positive" + TB + "mismatches" + TB
  + "raw-score" + TB + "bit-score" + TB + "e-score" + StringUtils.LS

  + "templateName" + TB + "+1" + TB + "0" + TB
  + "1" + TB + "12" + TB + "22" + TB + "1" + TB + "36" + TB + "36" + TB
  // matches   +matches   mismatches
  + "12" + TB + "100" + TB + "12" + TB + "100" + TB + "0" + TB
  // alignment-    bit-   e-scores
  + "-67" + TB + "30.4" + TB + "1.5e-8" + StringUtils.LS

  + "templateName" + TB + "+3" + TB + "1" + TB
  + "2" + TB + "12" + TB + "22" + TB + "3" + TB + "35" + TB + "36" + TB
  // matches   +matches   mismatches
  + "11" + TB + "100" + TB + "11" + TB + "100" + TB + "0" + TB
  + "-62" + TB + "28.5" + TB + "5.8e-8" +  StringUtils.LS;

  /**
   * Test method for {@link com.rtg.protein.ProteinOutputProcessor#process(long, java.lang.String, int, int, int, int)}.
   * @throws IOException if error
   * @throws InvalidParamsException
   */
  public final void testProcessNoProteinOutput() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tmp = FileUtils.createTempDir("proteinop", "filter");
    try {

      final File input = new File(tmp, "1");
      final SequencesReader r = ReaderTestUtils.getReaderDNA(READS_FASTA_PERFECT, input, null);
      r.close();
      final File input2 = new File(tmp, "2");
      final SequencesReader t = ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, input2);
      t.close();

      final NgsFilterParams filter = NgsFilterParams.builder()
      .outputFilter(OutputFilter.PROTEIN_TOPN)
      .useids(true)
      .errorLimit(CommonFlags.MAX_SCORE)
      .preFilterMinScore(0)
      .preFilterMinOverlap(0)
      .create();
      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).outputUnmapped(true).filterParams(filter).outputProteinSequences(false).create())
        .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).loadNames(true).create())
        .searchParams(SequenceParams.builder().directory(input2).useMemReader(true).loadNames(true).create())
        .proteinScoringMatrix(new ProteinScoringMatrix())
        .create();
      final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null);
      try {
        assertTrue(p instanceof ProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
      } finally {
        p.finish();
        p.close();
      }
      final String res = FileUtils.fileToString(new File(tmp, "alignments.tsv"));
      assertTrue(res.contains("#Version" + "\t" + Environment.getVersion()));
      assertTrue(res.contains(ProteinOutputProcessor.MAPX_OUTPUT_VERSION_HEADER + "\t" + ProteinOutputProcessor.MAPX_OUTPUT_VERSION));
      assertTrue(res.contains(ProteinOutputProcessor.MAPX_TEMPLATE_SDF_ID_HEADER + "\t" + t.getSdfId()));
      assertTrue(res.contains(ProteinOutputProcessor.MAPX_READ_SDF_ID_HEADER + "\t" + r.getSdfId()));
      MapXCliTest.checkAlignmentsNoHeader(EXPECTED_NO_PROTEIN, new File(tmp, "alignments.tsv"));
      MapXCliTest.checkUnmappedNoHeader(EXPECTED_UNMAPPED, new File(tmp, ProteinOutputProcessor.UNMAPPED_FILE));
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
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
