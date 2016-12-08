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

import java.io.File;
import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class TopEqualProteinOutputProcessorTest extends TestCase {

  /**
   */
  public TopEqualProteinOutputProcessorTest(String name) {
    super(name);
  }

  public final void testTopEqualProteinOutputProcessor() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tmp = FileUtils.createTempDir("proteinop", "filter");
    try {
    final NgsParams params = createParams(tmp, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 5);
    final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
    try {
     assertTrue(p instanceof TopEqualProteinOutputProcessor);
     p.process(0, "+1", 0, 0, 0, 0);
     p.process(0, "+1", 8, 0, 0, 0);
    } finally {
      p.finish();
      p.close();
    }

    MapXCliTest.checkAlignmentsNoHeader(ProteinOutputProcessorTest.EXPECTED_1, new File(tmp, "alignments.tsv"));
    MapXCliTest.checkUnmappedNoHeader(ProteinOutputProcessorTest.EXPECTED_UNMAPPED, new File(tmp, ProteinOutputProcessor.UNMAPPED_FILE));
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public static NgsParams createParams(final File tmp, String dna, String protein, final int topn) throws IOException, InvalidParamsException {
    final File input = new File(tmp, "1");
    ReaderTestUtils.getReaderDNA(dna, input, null).close();
    final File input2 = new File(tmp, "2");
    ReaderTestUtils.getReaderProtein(protein, input2).close();

    final NgsFilterParams filterParams = NgsFilterParams.builder().topN(topn).create();
    return NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).outputUnmapped(true).filterParams(filterParams).create())
                                                .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).loadNames(true).create())
                                                .searchParams(SequenceParams.builder().directory(input2).useMemReader(true).loadNames(true).create())
                                                .proteinScoringMatrix(new ProteinScoringMatrix())
                                                .create();
  }

  private static final String EXPECTED_2 =   "#template-name\tframe\tread-id\ttemplate-start\ttemplate-end\ttemplate-length\tread-start\tread-end\tread-length\ttemplate-protein\tread-protein\talignment\tidentical\t%identical\tpositive\t%positive\tmismatches\traw-score\tbit-score\te-score" + StringUtils.LS
    + "templateName\t+1\t0\t1\t12\t22\t1\t36\t36\tkwrknrkskknq\tkwrknrkskknq\tkwrknrkskknq\t12\t100\t12\t100\t0\t-67\t30.4\t1.5e-8" + StringUtils.LS
    + "templateName\t+3\t1\t2\t12\t22\t3\t35\t36\twrknrkskknq\twrknrkskknq\twrknrkskknq\t11\t100\t11\t100\t0\t-62\t28.5\t5.8e-8" + StringUtils.LS;

  public final void testTopEqualProteinOutputProcessor2() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tmp = FileUtils.createTempDir("proteinop", "filter");
    try {
      final NgsParams params = createParams(tmp, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 2);
      final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
      try {
        assertTrue(p instanceof TopEqualProteinOutputProcessor);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 0, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
      } finally {
        p.finish();
        p.close();
      }
      MapXCliTest.checkAlignmentsNoHeader(EXPECTED_2, new File(tmp, "alignments.tsv"));
      MapXCliTest.checkUnmappedNoHeader(ProteinOutputProcessorTest.EXPECTED_UNMAPPED, new File(tmp, ProteinOutputProcessor.UNMAPPED_FILE));
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public final void testCountOverflow() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tmp = FileUtils.createTempDir("proteinop", "filter");
    try {
      final NgsParams params = createParams(tmp, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 5);
      final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
      try {
        assertTrue(p instanceof TopEqualProteinOutputProcessor);
        for (int i = 0; i < 257; ++i) {
          p.process(0, "+1", 0, 0, 0, 0);
        }
        p.process(0, "+1", 8, 0, 0, 0);
        p.process(0, "+1", 8, 0, 0, 0);
      } finally {
        p.finish();
        p.close();
      }
      MapXCliTest.checkAlignmentsNoHeader(EXPECTED_2, new File(tmp, "alignments.tsv"));
      MapXCliTest.checkUnmappedNoHeader(ProteinOutputProcessorTest.EXPECTED_UNMAPPED, new File(tmp, ProteinOutputProcessor.UNMAPPED_FILE));
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testnvalidTopN() throws IOException, InvalidParamsException  {
    Diagnostic.setLogStream();
    final File tmp = FileUtils.createTempDir("proteinop", "invalidtopn");
    try {
      try {
      final NgsParams params = createParams(tmp, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 0);
      OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
      fail();
      } catch (final IllegalArgumentException ex) {
        //okay
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
    final File tmp2 = FileUtils.createTempDir("proteinop", "invalidtopn");
    try {
      try {
        final NgsParams params = createParams(tmp2, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 251);
        OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
        fail();
        } catch (final IllegalArgumentException ex) {
          //okay
        }
    } finally {
      assertTrue(FileHelper.deleteAll(tmp2));
    }
  }
}
