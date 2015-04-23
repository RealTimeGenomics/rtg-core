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
package com.rtg.variant.eval;

import java.io.File;
import java.util.List;

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class VcfRecordTabixCallableTest extends TestCase {

    public void testSomeMethod() throws Exception {
    final File dir = FileUtils.createTempDir("tabixVariantRunnerTest", "test");
    try {
      final File input = new File(dir, "snp_only.vcf.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz", input);
      final File tabix = new File(dir, "snp_only.vcf.gz.tbi");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz.tbi", tabix);
      ReferenceRanges ranges = new ReferenceRanges(true);
      final VcfRecordTabixCallable runner = new VcfRecordTabixCallable(input, ranges.forSequence("simulatedSequence13"), "simulatedSequence13", -1, VariantSetType.BASELINE, null, RocSortValueExtractor.NULL_EXTRACTOR, true, false, 100);
      List<DetectedVariant> set = runner.call().mVariants;
      assertEquals(2, set.size());
      final VcfRecordTabixCallable runner2 = new VcfRecordTabixCallable(input, ranges.forSequence("simulatedSequence2"), "simulatedSequence2", -1, VariantSetType.BASELINE, null, RocSortValueExtractor.NULL_EXTRACTOR, true, false, 100);
      set = runner2.call().mVariants;
      assertEquals(4, set.size());
      assertEquals(215, set.get(0).getStart());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
  public void testMissingSample() throws Exception {
    final MemoryPrintStream mp = new MemoryPrintStream();
    Diagnostic.setLogStream(mp.printStream());
    try {
      final File dir = FileUtils.createTempDir("tabixVariantRunnerTest", "test");
      try {
        final File input = new File(dir, "foo.vcf.gz");
        FileHelper.resourceToFile("com/rtg/sam/resources/vcf.txt.gz", input);
        final File tabix = new File(dir, "foo.vcf.gz.tbi");
        FileHelper.resourceToFile("com/rtg/sam/resources/vcf.txt.gz.tbi", tabix);
        ReferenceRanges ranges = new ReferenceRanges(true);
        final VcfRecordTabixCallable runner = new VcfRecordTabixCallable(input, ranges.forSequence("20"), "20", -1, VariantSetType.CALLS, "asdf", RocSortValueExtractor.NULL_EXTRACTOR, true, false, 100);
        try {
          runner.call();
          fail();
        } catch (NoTalkbackSlimException e) {
          TestUtils.containsAll(e.toString(), "Sample \"asdf\" not found in calls VCF");
        }

      } finally {
        FileHelper.deleteAll(dir);
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

}
