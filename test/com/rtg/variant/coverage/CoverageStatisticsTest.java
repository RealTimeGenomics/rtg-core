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

  public void testFold80() {
    final long[] hist = {6293, 1281, 2007, 869, 1433, 811, 1177, 699, 956, 620, 800, 906, 763, 664, 907, 685, 776, 665, 768, 790, 690, 750, 730, 726, 640, 661, 837, 758, 643, 595, 584, 500, 582, 579, 700, 706, 622, 611, 614, 615, 591, 672, 615, 583, 583, 478, 582, 449, 613, 582, 618, 504, 504, 591, 565, 670, 587, 465, 619, 631, 650, 672, 646, 543, 584, 605, 470, 584, 556, 610, 769, 599, 682, 722, 559, 587, 587, 539, 782, 735, 546, 631, 419, 555, 558, 651, 664, 706, 653, 562, 681, 592, 607, 603, 655, 644, 598, 607, 652, 705, 754, 682, 724, 695, 589, 627, 651, 640, 569, 643, 753, 734, 706, 543, 621, 655, 606, 557, 557, 478, 590, 556, 562, 595, 617, 672, 652, 611, 662, 646, 719, 799, 752, 628, 709, 737, 557, 575, 650, 684, 654, 718, 692, 735, 823, 751, 822, 742, 879, 768, 844, 720, 798, 769, 676, 770, 780, 769, 737, 717, 761, 779, 827, 900, 846, 921, 843, 822, 795, 874, 912, 943, 847, 854, 817, 864, 837, 893, 923, 994, 978, 873, 903, 888, 914, 839, 843, 783, 920, 832, 880, 859, 934, 848, 849, 839, 919, 980, 903, 970, 873, 847, 948, 1031, 993, 1000, 938, 912, 924, 991, 900, 911, 958, 907, 994, 972, 955, 994, 924, 1016, 974, 959, 969, 1098, 1128, 1202, 1113, 1127, 1050, 1148, 1195, 1164, 1094, 1084, 1153, 1008, 1124, 1107, 1212, 1165, 1170, 1202, 1487, 1398, 1365, 1299, 1210, 1238, 1239, 1205, 1185, 1165, 1261, 1196, 1319, 1298, 1338, 1334, 1171, 1089, 1199, 1239, 1285, 1156, 1247, 1251, 1273, 1244, 1195, 1234, 1222, 1329, 1283, 1354, 1369, 1314, 1377, 1309, 1259, 1225, 1337, 1274, 1432, 1373, 1375, 1337, 1351, 1305, 1332, 1233, 1312, 1336, 1320, 1310, 1458, 1387, 1350, 1390, 1318, 1282, 1327, 1382, 1194, 1204, 1311, 1305, 1307, 1238, 1258, 1253, 1308, 1298, 1321, 1239, 1287, 1359, 1291, 1302, 1284, 1255, 1246, 1222, 1175, 1157, 1142, 1257, 1200, 1214, 1289, 1193, 1267, 1185, 1268, 1237, 1206, 1231, 1095, 1163, 1160, 1153, 1209, 1139, 1107, 1093, 1110, 1058, 1131, 1205, 1200, 1068, 1078, 1022, 1025, 1016, 966, 927, 953, 928, 912, 814, 942, 925, 896, 910, 859, 954, 882, 843, 928, 873, 817, 871, 882, 867, 822, 863, 826, 801, 828, 704, 737, 759, 750, 809, 762, 706, 730, 768, 772, 736, 755, 678, 701, 715, 708, 664, 679, 642, 689, 631, 629, 597, 566, 617, 596, 584, 545, 541, 579, 546, 565, 539, 538, 542, 529, 548, 502, 554, 540, 536, 555, 512, 487, 456, 471, 452, 484, 508, 458, 482, 560, 544, 524, 465, 423, 452, 441, 439, 458, 419, 402, 394, 426, 441, 378, 421, 402, 382, 437, 399, 364, 386, 328, 328, 321, 334, 289, 330, 289, 290, 277, 332, 360, 311, 322, 265, 259, 255, 270, 246, 231, 246, 236, 225, 240, 248, 267, 271, 233, 201, 268, 192, 264, 202, 256, 204, 187, 210, 214, 200, 206, 172, 188, 186, 202, 213, 182, 213, 182, 141, 196, 171, 186, 208, 263, 238, 226, 221, 190, 150, 149, 142, 117, 109, 112, 124, 138, 154, 140, 147, 154, 158, 154, 134, 149, 137, 126, 110, 134, 138, 115, 87, 98, 101, 84, 106, 103, 105, 106, 124, 121, 92, 131, 106, 86, 104, 83, 98, 81, 82, 68, 78, 87, 64, 75, 75, 71, 70, 50, 67, 62, 63, 60, 58, 69, 54, 46, 48, 48, 65, 81, 53, 56, 49, 55, 32, 44, 55, 68, 65, 64, 61, 56, 45, 37, 53, 43, 29, 45, 51, 53, 31, 44, 44, 39, 37, 48, 41, 50, 50, 36, 50, 30, 40, 27, 29, 28, 45, 42, 44, 26, 30, 32, 32, 30, 30, 15, 35, 40, 31, 46, 39, 39, 25, 26, 22, 11, 23, 24, 30, 32, 27, 49};

    final Double f = CoverageStatistics.fold80(419499, 419499 * 245.54, hist);
    assertNotNull(f);
    assertEquals(2.11, f, 0.01);

    final long[] hist2 = {6293};
    final Double f2 = CoverageStatistics.fold80(419499, 419499 * 5010.8, hist2);
    assertNull(f2);
  }

  public void testFold80b() {
    final long[] hist3 = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    assertEquals(4.0, CoverageStatistics.fold80(100, 100 * 8, hist3), 0.01);
    hist3[1] = 9;
    assertEquals(3.8, CoverageStatistics.fold80(100, 100 * 8, hist3), 0.01);
    hist3[1] = 11;
    assertEquals(4.2, CoverageStatistics.fold80(100, 100 * 8, hist3), 0.01);
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
