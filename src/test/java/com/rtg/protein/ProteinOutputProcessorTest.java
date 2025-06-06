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
import com.rtg.launcher.ISequenceParams;
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

  // AAA TGG CGC AAA AAC AGA AAG TCG AAA AAA AAT CAA AGA AAT TAT AAC CAC GAC GCA GCA GAC GCA
  // 0   1   2   3   4   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20
  private static final String TEMPLATE_DNA = "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGAAATTATAACCACGACGCAGCAGACGCA";
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

  private static final String[] PROTEIN_QUERY_PERFECT = {
    TestUtils.dnaToProtein("AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA"),
    TestUtils.dnaToProtein("TGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGA"),
    TestUtils.dnaToProtein("AATCAAAGAAATTATAACCACGACGCAGCAGACGCA"),
  };
  static final String PROTEIN_QUERY_FASTA_PERFECT = toFasta("testProt", PROTEIN_QUERY_PERFECT);

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

  private static NgsParams createParams(File tmp, boolean outputprot) throws IOException {
    return createParams(tmp, outputprot, getSequenceParamsDna(READS_FASTA_PERFECT));
  }

  private static NgsParams createParams(File tmp, boolean outputprot, ISequenceParams buildparams) throws IOException {
    final NgsFilterParams filter = NgsFilterParams.builder()
    .outputFilter(OutputFilter.PROTEIN_TOPN)
    .useids(true)
    .errorLimit(MapFlags.MAX_SCORE)
    .preFilterMinScore(0)
    .preFilterMinOverlap(0)
    .create();
    return NgsParams.builder()
      .outputParams(NgsOutputParams.builder().outputDir(tmp).outputUnmapped(true).filterParams(filter).outputProteinSequences(outputprot).create())
      .buildFirstParams(buildparams)
      .searchParams(getSequenceParamsProtein(TEMPLATE_FASTA))
      .proteinScoringMatrix(new ProteinScoringMatrix())
      .create();
  }

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

  public final void testProcessProteinBuild() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteinop")) {
      final NgsParams params = createParams(tmp, true, getSequenceParamsProtein(PROTEIN_QUERY_FASTA_PERFECT));
      try (final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null)) {
        assertTrue(p instanceof ProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 1, 0, 0, 0);
        p.process(0, "+1", 2, 9, 0, 0);
        p.finish();
      }
      mNano.check("pop-prot-alignments.tsv", getFileFromMatch(new File(tmp, TABULAR_ALIGNMENTS), "#template-name"));
      mNano.check("pop-prot-unmapped.tsv", getFileFromMatch(new File(tmp, UNMAPPED_FILE), "#read-id"));
    }
  }

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
