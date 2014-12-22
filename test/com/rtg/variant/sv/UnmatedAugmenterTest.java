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

package com.rtg.variant.sv;

import java.io.File;
import java.io.InputStream;

import com.rtg.util.Resources;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class UnmatedAugmenterTest extends TestCase {

  private File mFile;
  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mFile = FileUtils.createTempDir("test", "t");
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mFile));
    mFile = null;
  }


  public void testMain() throws Exception {
    final File in = new File(mFile, "in.sam");
    final File out = new File(mFile, "out.sam");
    InputStream stream = Resources.getResourceAsStream("com/rtg/sam/resources/unmated.sam");
    try {
      FileHelper.streamToFile(stream, in);
    } finally {
      stream.close();
    }
    final String exp;
    stream = Resources.getResourceAsStream("com/rtg/sam/resources/augmented.sam");
    try {
      exp = FileUtils.streamToString(stream);
    } finally {
      stream.close();
    }
    final UnmatedAugmenter un = new UnmatedAugmenter();
    un.augmentUnmated(in, out, false, null);
    final String outStr = FileUtils.fileToString(out);
    final String expNoPg = exp.replaceAll("@PG.*\n", "");
    final String outStrNoPg = outStr.replaceAll("@PG.*\n", "");
    assertEquals(expNoPg, outStrNoPg);
    assertTrue(outStr, outStr.contains("@PG\tID:rtg"));
  }

  public void testAugmenting() throws Exception {
    final File mated = new File(mFile, "mated.sam");
    final File unmated = new File(mFile, "unmated.sam");
    final File unmapped = new File(mFile, "unmapped.sam");
    final File outunmated = new File(mFile, "outunmated.sam");
    final File outunmapped = new File(mFile, "outunmapped.sam");


    InputStream stream = Resources.getResourceAsStream("com/rtg/sam/resources/mergemated.sam");
    try {
      FileHelper.streamToFile(stream, mated);
    } finally {
      stream.close();
    }
    stream = Resources.getResourceAsStream("com/rtg/sam/resources/mergeunmated.sam");
    try {
      FileHelper.streamToFile(stream, unmated);
    } finally {
      stream.close();
    }
    stream = Resources.getResourceAsStream("com/rtg/sam/resources/mergeunmapped.sam.gz");
    try {
      FileHelper.streamToFile(stream, unmapped);
    } finally {
      stream.close();
    }
    final String exp;
    stream = Resources.getResourceAsStream("com/rtg/sam/resources/mergeunmapped-aug.sam");
    try {
      exp = FileUtils.streamToString(stream);
    } finally {
      stream.close();
    }
    final UnmatedAugmenter un = new UnmatedAugmenter();
    final ReadGroupStatsCalculator calc = new ReadGroupStatsCalculator();
    UnmatedAugmenter.populateReadGroupStats(mated, calc);
    un.augmentUnmated(unmated, outunmated, false, calc);
    un.augmentUnmapped(unmapped, outunmapped, false, calc);
    final String outUnmappedStr = FileUtils.fileToString(outunmapped);
    final String expNoPg = exp.replaceAll("@PG.*\n", "");
    final String outStrNoPg = outUnmappedStr.replaceAll("@PG.*\n", "");
    assertEquals(expNoPg, outStrNoPg);
    assertTrue(outUnmappedStr, outUnmappedStr.contains("@PG\tID:slim"));
  }

  private SAMRecord createSAMRecord(SAMFileHeader header, String readName, int flags, int position) {
    final SAMRecord rec = new SAMRecord(header);
    rec.setAlignmentStart(position);
    rec.setFlags(flags);
    rec.setReadName(readName);
    rec.setReferenceName("simulatedSequence10");
    rec.setCigarString("20=");
    rec.setReadString("AAAAAAAAAAAAAAAAAAAA");
    rec.setBaseQualityString("CGHFDCI7GHIIDH9?I?CC");
    rec.setMappingQuality(37);
    rec.setAttribute("AS", 0);
    rec.setAttribute("NH", 1);
    rec.setAttribute("IH", 1);
    rec.setAttribute("RG", "rg1");
    rec.setMateReferenceName("*");
    rec.setMateAlignmentStart(0);
    rec.setInferredInsertSize(0);
    return rec;
  }

  public void testBlend() {
    final SAMFileHeader header = new SAMFileHeader();
    header.addSequence(new SAMSequenceRecord("simulatedSequence10", 60977));
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
    rg.setPlatform("ILLUMINA");
    rg.setSample("sm1");
    header.addReadGroup(rg);
    final SAMRecord sr1a = createSAMRecord(header, "96784", 129, 5214);
    final SAMRecord sr1b = createSAMRecord(header, "96784", 81, 6128);
    final SAMRecord sr2a = createSAMRecord(header, "42112", 65, 5186);
    final SAMRecord sr2b = createSAMRecord(header, "42112", 145, 6125);
    final SAMRecord sr3 = createSAMRecord(header, "98275", 137, 5108);
    final UnmatedAugmenter.Merger merger = new UnmatedAugmenter.Merger();
    final UnmatedAugmenter ua1 = merger.createUnmatedAugmenter();
    ua1.processRecord(sr1a);
    ua1.processRecord(sr2b);
    final UnmatedAugmenter ua2 = merger.createUnmatedAugmenter();
    ua2.processRecord(sr2a);
    ua2.processRecord(sr1b);
    ua2.processRecord(sr1a);
    ua2.processRecord(sr3);

    final UnmatedAugmenter uaBlend = merger.blend();
    assertEquals(uaBlend, merger.blend());

    uaBlend.updateUnmatedRecord(sr1a);
    uaBlend.updateUnmatedRecord(sr1b);
    uaBlend.updateUnmatedRecord(sr2a);
    uaBlend.updateUnmatedRecord(sr2b);
    uaBlend.updateUnmatedRecord(sr3);
    checkRecord(sr3, 137, "*", 0, 0);
    checkRecord(sr1a, 161, "simulatedSequence10", 6128, 934);
    checkRecord(sr1b, 81, "simulatedSequence10", 5214, -934);
    checkRecord(sr2a, 97, "simulatedSequence10", 6125, 959);
    checkRecord(sr2b, 145, "simulatedSequence10", 5186, -959);
  }

  private void checkRecord(SAMRecord rec, int flags, String mateReferenceName, int alignmentStart, int insertSize) {
    assertEquals(flags, rec.getFlags());
    assertEquals(mateReferenceName, rec.getMateReferenceName());
    assertEquals(alignmentStart, rec.getMateAlignmentStart());
    assertEquals(insertSize, rec.getInferredInsertSize());
  }

  public void testGz() throws Exception {
    final File in = new File(mFile, "in.sam.gz");
    final File out = new File(mFile, "out.sam.gz");
    FileHelper.resourceToGzFile("com/rtg/sam/resources/unmated.sam", in);
    final String exp;
    exp = FileHelper.resourceToString("com/rtg/sam/resources/augmented.sam");
    final UnmatedAugmenter un = new UnmatedAugmenter();
    un.augmentUnmated(in, out, true, null);
    final String outStr = FileHelper.gzFileToString(out);
    final String expNoPg = exp.replaceAll("@PG.*\n", "");
    final String outStrNoPg = outStr.replaceAll("@PG.*\n", "");
    assertEquals(expNoPg, outStrNoPg);
    assertTrue(outStr.contains("@PG\tID:rtg"));
  }
}
