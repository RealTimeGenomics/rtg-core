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

import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsParams;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class TopNProteinOutputProcessorTest extends TestCase {

  /**
   */
  public TopNProteinOutputProcessorTest(String name) {
    super(name);
  }

 private static final String EXPECTED = ""
   + "#template-name\tframe\tread-id\ttemplate-start\ttemplate-end\ttemplate-length\t"
   + "read-start\tread-end\tread-length\ttemplate-protein\tread-protein\talignment\tidentical\t%identical\tpositive\t%positive\tmismatches\traw-score\tbit-score\te-score" + StringUtils.LS
   + null + "\t+1\t0\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t1\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t1\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t2\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t2\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t3\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t3\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t4\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t4\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t4\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t6\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS
   + null + "\t+1\t6\t1\t0\t-1\t1\t0\t-1\t\t\t\t0\t0\t0\t0\t0\t0\t4.6\t0" + StringUtils.LS;

  public void testTopN() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tmp2 = FileUtils.createTempDir("proteintopn", "filter");
    try {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp2, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 3);
      final TopNProteinOutputProcessor p = new TopNProteinOutputProcessor(params, null);
      try {
        p.writeResult(getProteinRes(0, -42));
        p.writeResult(getProteinRes(1, 42));
        p.writeResult(getProteinRes(1, -42));
        p.writeResult(getProteinRes(2, -42));
        p.writeResult(getProteinRes(2, 42));
        p.writeResult(getProteinRes(3, -42));
        p.writeResult(getProteinRes(3, -42));
        p.writeResult(getProteinRes(4, -42));
        p.writeResult(getProteinRes(4, -42));
        p.writeResult(getProteinRes(4, -42));
        p.writeResult(getProteinRes(5, -42));
        p.writeResult(getProteinRes(5, -42));
        p.writeResult(getProteinRes(5, -42));
        p.writeResult(getProteinRes(5, -42));
        p.writeResult(getProteinRes(6, 9));
        p.writeResult(getProteinRes(6, 2));
        p.writeResult(getProteinRes(6, 2));
        p.writeResult(getProteinRes(6, 3));
        p.writeResult(getProteinRes(6, 3));
      } finally {
        p.finish();
        p.close();
      }
      MapXCliTest.checkAlignmentsNoHeader(EXPECTED, new File(tmp2, "alignments.tsv"));
    } finally {
      assertTrue(FileHelper.deleteAll(tmp2));
    }
  }

  public void testTopN250() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File tmp2 = FileUtils.createTempDir("proteintopn", "filter");
    try {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp2, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 250);
      final TopNProteinOutputProcessor p = new TopNProteinOutputProcessor(params, null);
      try {
        for (int k = 0; k < 250; k++) {
          p.writeResult(getProteinRes(0, -47));
        }
        for (int k = 0; k < 251; k++) {
          p.writeResult(getProteinRes(1, 47));
        }
      } finally {
        p.finish();
        p.close();
      }
      MapXCliTest.checkAlignmentsNoHeader(StringUtils.convertLineEndings(FileHelper.resourceToString("com/rtg/protein/resources/250.txt")), new File(tmp2, "alignments.tsv"));
    } finally {
      assertTrue(FileHelper.deleteAll(tmp2));
    }
  }

  private ProteinAlignmentResult getProteinRes(int readid, final int score) throws InvalidParamsException, IOException {
    final SharedProteinResources resx = new SharedProteinResources(new ProteinScoringMatrix(), null, null, false);
    return new ProteinAlignmentResult(resx, 0, readid * 6, (int[]) null, 0) {
        @Override
        int alignmentScore() {
          return score;
        }
      };
  }

}
