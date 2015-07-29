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
package com.rtg.variant.coverage;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.GlobalFlags;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.IndexUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Date: 28/11/12
 * Time: 5:06 PM
 */
public class CoverageStatisticsTest extends TestCase {

  public void testAddAverageCoverage() {
    final CoverageStatistics cs = new CoverageStatistics(null, false);

    final RangeList.RangeData<String> range = new RangeList.RangeData<>(0, 10, "blah");
    final List<RangeList.RangeData<String>> ranges = new ArrayList<>();
    ranges.add(range);
    final RangeList<String> rangesearch = new RangeList<>(ranges);
    cs.setRange(rangesearch.getRangeList().get(0));

//    cs.addCoveredBasesCount("blah", 5);
//    cs.addNonNCoveredBasesCount("blah", 3);
//
//    cs.addAverageCoverage("blah", 5, 3, 10, 8);
    cs.updateCoverageHistogram(1, false, 1);
    cs.updateCoverageHistogram(1, false, 1);
    cs.updateCoverageHistogram(1, false, 1);
    cs.updateCoverageHistogram(1, true, 1);
    cs.updateCoverageHistogram(1, true, 1);
    for (int i = 0; i < 5; i++) {
      cs.updateCoverageHistogram(0, false, 1);
    }

    cs.setRange(null);
    final String exp = "Coverage per region:" + LS
                         + "   depth  breadth  covered  size         name" + LS
                         + "  0.3750   0.3750        3     8         blah" + LS
                         + "  0.3750   0.3750        3     8  all regions" + LS;
    assertEquals(exp, cs.getStatistics());
  }

  public void testBedRegions() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File samFile = new File(tmpDir, "sam.sam.gz");
      IndexUtils.ensureBlockCompressed(FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", samFile));
      new TabixIndexer(samFile).saveSamIndex();

      final String bedRegion = ""
                                       + "simulatedSequence1\t40\t47\n"
                                       + "simulatedSequence1\t55\t65\tblah\n"
                                       + "simulatedSequence1\t60\t75\tfeh\n"
                                       + "simulatedSequence2\t17\t43\n"
                                       + "simulatedSequence2\t55\t59\tfeh\n"
                                       + "simulatedSequence2\t63\t94\n"
              ;
      final File bedRegionsFile = FileHelper.stringToGzFile(bedRegion, new File(tmpDir, "bedRegions.bed.gz"));

      final File output = new File(tmpDir, "output");

      try (MemoryPrintStream mps = new MemoryPrintStream()) {
        final String tmpl = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAA";

        final File template = ReaderTestUtils.getDNADir(">simulatedSequence1\n" + tmpl + "\n>simulatedSequence2\n" + tmpl + "\n", new File(tmpDir, "template"));
        GlobalFlags.resetAccessedStatus();
        final int code = new CoverageCli().mainInit(new String[]{"-o", output.getPath(), "-s", "0", samFile.getPath(),
                                                                        "--bed-regions", bedRegionsFile.getPath(), "-t", template.getPath()}, mps.outputStream(), mps.printStream());
        assertEquals(mps.toString(), 0, code);

        final String summary = FileHelper.fileToString(new File(output, "summary.txt"));
        final String expSumm = "Coverage per region:\n"
                                       + "   depth  breadth  covered  size                      name\n"
                                       + "  2.0000   1.0000        7     7  simulatedSequence1:41-47\n"
                                       + "  2.1111   1.0000        9     9                      blah\n"
                                       + "  2.4118   1.0000       17    17                       feh\n"
                                       + "  1.7692   1.0000       26    26  simulatedSequence2:18-43\n"
                                       + "  1.8000   0.7333       22    30  simulatedSequence2:64-94\n"
                                       + "  1.9294   0.9059       77    85               all regions\n";
        assertEquals(expSumm.replaceAll("\n|\r\n", LS), summary);

        final String stats = FileHelper.fileToString(new File(output, "stats.tsv"));
        final String expStats = "#depth\tbreadth\tcovered\tsize\tname\n"
                                        + "2.0000\t1.0000\t7\t7\tsimulatedSequence1:41-47\n"
                                        + "2.1111\t1.0000\t9\t9\tblah\n"
                                        + "2.4118\t1.0000\t17\t17\tfeh\n"
                                        + "1.7692\t1.0000\t26\t26\tsimulatedSequence2:18-43\n"
                                        + "1.8000\t0.7333\t22\t30\tsimulatedSequence2:64-94\n"
                                        + "1.9294\t0.9059\t77\t85\tall regions\n";
        assertEquals(expStats.replaceAll("\n|\r\n", LS), stats);

        final String levels = FileHelper.fileToString(new File(output, "levels.tsv"));
        final String explevels = "#coverage_level\tcount\t%age\t%cumulative\n"
                                         + "0\t8\t9.41\t100.00\n"
                                         + "1\t13\t15.29\t90.59\n"
                                         + "2\t46\t54.12\t75.29\n"
                                         + "3\t13\t15.29\t21.18\n"
                                         + "4\t5\t5.88\t5.88\n";
        assertEquals(explevels.replaceAll("\n|\r\n", LS), levels);

        final String report = FileHelper.fileToString(new File(output, "index.html"));
        assertTrue(report.contains("0</td><td align=\"right\">9.41</td>"));
      }
    }
  }
}
