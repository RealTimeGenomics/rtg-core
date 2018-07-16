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
import java.io.IOException;
import java.io.StringReader;
import java.util.HashMap;
import java.util.Map;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.Sex;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.PairOrientation;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Test class
 */
public class UnmatedAugmenterTest extends AbstractNanoTest {

  public void testMain() throws Exception {
    try (TestDirectory temp = new TestDirectory()) {
      final File in = new File(temp, "in.sam");
      final File out = new File(temp, "out.sam");
      FileHelper.resourceToFile("com/rtg/sam/resources/unmated.sam", in);
      final UnmatedAugmenter un = new UnmatedAugmenter();
      un.augmentUnmated(in, out, new ReadGroupStatsCalculator());
      final String outStr = FileUtils.fileToString(out);
      final String outStrNoPg = StringUtils.grepMinusV(outStr, "^@PG");
      mNano.check("augmented.sam", outStrNoPg);
    }
  }

  public void testGz() throws Exception {
    try (TestDirectory temp = new TestDirectory()) {
      final File in = new File(temp, "in.sam.gz");
      final File out = new File(temp, "out.sam.gz");
      FileHelper.resourceToGzFile("com/rtg/sam/resources/unmated.sam", in);
      final UnmatedAugmenter un = new UnmatedAugmenter();
      un.augmentUnmated(in, out, new ReadGroupStatsCalculator());
      final String outStr = FileHelper.gzFileToString(out);
      final String outStrNoPg = outStr.replaceAll("@PG.*\n", "");
      mNano.check("augmented.sam", outStrNoPg);
    }
  }

  public void testAugmenting() throws Exception {
    try (TestDirectory temp = new TestDirectory()) {
      final File mated = new File(temp, "mated.sam");
      final File unmated = new File(temp, "unmated.sam");
      final File unmapped = new File(temp, "unmapped.sam");
      final File outunmated = new File(temp, "outunmated.sam");
      final File outunmapped = new File(temp, "outunmapped.sam");

      FileHelper.resourceToFile("com/rtg/sam/resources/mergemated.sam", mated);
      FileHelper.resourceToFile("com/rtg/sam/resources/mergeunmated.sam", unmated);
      FileHelper.resourceToFile("com/rtg/sam/resources/mergeunmapped.sam.gz", unmapped);


      final UnmatedAugmenter un = new UnmatedAugmenter();
      final ReadGroupStatsCalculator calc = new ReadGroupStatsCalculator();
      calc.addFile(mated);
      un.augmentUnmated(unmated, outunmated, calc);
      un.augmentUnmapped(unmapped, outunmapped, calc);
      final String outUnmappedStr = FileUtils.fileToString(outunmapped);
      final String outStrNoPg = StringUtils.grepMinusV(outUnmappedStr, "^@PG");
      mNano.check("mergeunmapped-aug.sam", outStrNoPg);
    }
  }

  public void testAugmentingMixed() throws Exception {
    try (TestDirectory temp = new TestDirectory()) {
      final File alignments = new File(temp, "alignments.sam.gz");
      final File outalignments = new File(temp, "outalignments.sam.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/mergecombined.sam.gz", alignments);
      final ReadGroupStatsCalculator calc = new ReadGroupStatsCalculator();
      new UnmatedAugmenter().augmentMixed(alignments, outalignments, calc);
      final String outUnmappedStr = FileHelper.gzFileToString(outalignments);
      final String outStrNoPg = StringUtils.grepMinusV(outUnmappedStr, "^@PG");
      mNano.check("mergecombined-aug.sam", outStrNoPg);
    }
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
    ua1.addRecord(sr1a);
    ua1.addRecord(sr2b);
    final UnmatedAugmenter ua2 = merger.createUnmatedAugmenter();
    ua2.addRecord(sr2a);
    ua2.addRecord(sr1b);
    ua2.addRecord(sr1a);
    ua2.addRecord(sr3);

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


  private static final String PSEUDO_REF_2 = "version 1\n"
    + "either\tdef\tdiploid\tlinear\n"
    + "male    seq     chrX    haploid linear  chrY\n"
    + "male    seq     chrY    haploid linear  chrX\n"
    + "female  seq     chrX    diploid linear\n"
    + "female  seq     chrY    none    linear\n"
    + "#PAR1 pseudoautosomal region\n"
    + "male    dup     chrX:10001-2781479      chrY:10001-2781479\n"
    + "#PAR2 pseudoautosomal region\n"
    + "male    dup     chrX:155701383-156030895        chrY:56887903-57217415\n";

  public ReferenceGenome makeReferenceGenome() throws IOException {
    final Map<String, Integer> names = new HashMap<>();
    names.put("chrX", 156040895);
    names.put("chrY", 57227415);
    return new ReferenceGenome(names, new StringReader(PSEUDO_REF_2), Sex.MALE);
  }

  private SAMRecord makeRec(SAMFileHeader header) {
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadString(StringUtils.repeat("a", 150));
    rec.setReadUnmappedFlag(true);
    rec.setMateUnmappedFlag(false);
    rec.setReadName("0");
    rec.setReferenceName("chrY");
    return rec;
  }

  public void testUnmappedPlacement() throws IOException {
    final SAMFileHeader header = new SAMFileHeader();
    header.addSequence(new SAMSequenceRecord("chrX", 156040895));
    header.addSequence(new SAMSequenceRecord("chrY", 57227415));

    final ReferenceGenome r = makeReferenceGenome();

    int mateStart = 2781679;
    int mateEnd = 2781830;
    final int parEnd = 2781479; // One based inclusive end position of the PAR region
    for (int i = -1; i < 2; i++) {
      SAMRecord rec = makeRec(header);
      final int expStart = parEnd - i;
      final int flen = 350 + i;
      final String expChr = expStart <= parEnd ? "chrX" : "chrY";
      UnmatedAugmenter.setPlacement(rec, null, MachineOrientation.FR, PairOrientation.R1, mateStart, mateEnd, flen);
      assertEquals("chrY", rec.getReferenceName());
      assertEquals(expStart, rec.getAlignmentStart());
      rec = makeRec(header);
      UnmatedAugmenter.setPlacement(rec, r, MachineOrientation.FR, PairOrientation.R1, mateStart, mateEnd, flen);
      assertEquals(expChr, rec.getReferenceName());
      assertEquals(expStart, rec.getAlignmentStart());
    }

    mateStart = 9799;
    mateEnd = 9950;
    final int parStart = 10001; // One based start position of the PAR region
    for (int i = -1; i < 2; i++) {
      SAMRecord rec = makeRec(header);
      final int expStart = parStart + i;
      final int flen = 351 + i;
      final String expChr = expStart >= parStart ? "chrX" : "chrY";
      UnmatedAugmenter.setPlacement(rec, null, MachineOrientation.FR, PairOrientation.F1, mateStart, mateEnd, flen);
      assertEquals("chrY", rec.getReferenceName());
      assertEquals(expStart, rec.getAlignmentStart());
      rec = makeRec(header);
      UnmatedAugmenter.setPlacement(rec, r, MachineOrientation.FR, PairOrientation.F1, mateStart, mateEnd, flen);
      assertEquals(expChr, rec.getReferenceName());
      assertEquals(expStart, rec.getAlignmentStart());
    }
  }

}
