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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.IndexUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

public class CoverageStatisticsTest extends AbstractNanoTest {

  public void testAddAverageCoverage() throws IOException {
    final CoverageStatistics cs = new CoverageStatistics(null, false);

    final RangeList.RangeData<String> range = new RangeList.RangeData<>(0, 10, "blah");
    final List<RangeList.RangeData<String>> ranges = new ArrayList<>();
    ranges.add(range);
    final RangeList<String> rangesearch = new RangeList<>(ranges);
    cs.setRange("seq", rangesearch.getRangeList().get(0));

    cs.updateCoverageHistogram(10, false, 1);
    cs.updateCoverageHistogram(10, false, 1);
    cs.updateCoverageHistogram(10, false, 1);
    cs.updateCoverageHistogram(10, true, 1);
    cs.updateCoverageHistogram(10, true, 1);
    for (int i = 0; i < 5; ++i) {
      cs.updateCoverageHistogram(0, false, 1);
    }

    cs.setRange("seq", null);
    mNano.check("covstats-add-median.txt", cs.getStatistics());
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

      final String tmpl = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAA";

      final File template = ReaderTestUtils.getDNADir(">simulatedSequence1\n" + tmpl + "\n>simulatedSequence2\n" + tmpl + "\n", new File(tmpDir, "template"));
      final MainResult r = MainResult.run(new CoverageCli(), "-o", output.getPath(), "-s", "0", samFile.getPath(),
        "--bed-regions", bedRegionsFile.getPath(), "-t", template.getPath());
      assertEquals(r.err(), 0, r.rc());

      mNano.check("covstats-summary.txt", FileHelper.fileToString(new File(output, "summary.txt")));
      mNano.check("covstats-stats.tsv", FileHelper.fileToString(new File(output, "stats.tsv")));
      mNano.check("covstats-levels.tsv", FileHelper.fileToString(new File(output, "levels.tsv")));

      final String report = FileHelper.fileToString(new File(output, "index.html"));
      assertTrue(report.contains("0</td><td align=\"right\">9.41</td>"));
    }
  }
}
