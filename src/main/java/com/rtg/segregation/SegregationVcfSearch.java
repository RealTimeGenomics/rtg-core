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

package com.rtg.segregation;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.rtg.bed.BedWriter;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.Pair;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.GzipAsynchOutputStream;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.ContigField;
import com.rtg.vcf.header.SampleField;
import com.rtg.vcf.header.VcfHeader;

/**
 * Analyse the platinum cohort given two parents and some number of children (11 for the platinum data set) from a VCF file.
 * Outputs the phased grouped results and those calls considered errors.
 * <br><br>
 * The prefixes on the individual lines are:
 * <ul>
 * <LI>OK - no errors or other problems - original line. Contributes to a block</LI>
 * <LI>HO - both parents homozygous, no phasing information - original line</LI>
 * <LI>NO - both parents no ploidy, no phasing information - original line</LI>
 * <LI>AH - parents and children all heterozygous (with same genotype). Very unlikely and provides no phasing information - original line</LI>
 * <LI>ME - Mendelian error - original line</LI>
 * <LI>BN - new block (at start or when unable to bridge a gap) - block</LI>
 * <LI>BL - phased block - block</LI>
 * <LI>BE - unphased block - block</LI>
 * <LI>BX - block after a cross over - block</LI>
 * <LI>NS - Start of new sequence - sequence name</LI>
 * <LI>XO - cross over - <code>fa(ther)</code>/<code>mo(ther)</code>, child which changed, level before, level after</LI>
 * </ul>
 */
public class SegregationVcfSearch extends AbstractCli {

  private static final String MODULE_NAME = "phasingsearch";

  private static final String VCF_FLAG = "vcf";
  private static final String OUTPUT_FLAG = "output";
  private static final String OUTPUT_REGIONS_FLAG = "regions-output";
  private static final String FATHER_FLAG = "father";
  private static final String MOTHER_FLAG = "mother";
  private static final String NEW_PENALTY_FLAG = "Xnew-penalty";
  private static final String XO_PENALTY_FLAG = "Xxo-penalty";

  private PrintStream mOut = null;
  private BedWriter mBed = null;
  private Search mSearch = null;
  private SegregationBlock mBlock = null;
  private int mLastLength;
  private String mLastSeq = null;
  final Map<String, Integer> mChrLengths = new HashMap<>();
  private Map<Pair<Sex, String>, ReferenceSequence> mPloidyMap;
  int mNewPenalty;
  int mXoPenalty;

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "search for phasing based on segregation analysis";
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);

    mFlags.setDescription("Analyze a VCF to produce phasing segregation information. VCF must only contain samples for the parents and the children of those parents.");

    CommonFlags.initReferenceTemplate(mFlags, true);
    mFlags.registerRequired(VCF_FLAG, File.class, CommonFlags.FILE, "input VCF file to be analyzed").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(OUTPUT_FLAG, File.class, CommonFlags.FILE, "output text file containing block and crossover information").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(OUTPUT_REGIONS_FLAG, File.class, CommonFlags.FILE, "output BED file to contain regions of phasing segregation information").setCategory(INPUT_OUTPUT);

    mFlags.registerRequired(FATHER_FLAG, String.class, CommonFlags.STRING, "sample name of the father").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(MOTHER_FLAG, String.class, CommonFlags.STRING, "sample name of the mother").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(NEW_PENALTY_FLAG, Integer.class, CommonFlags.INT, "override penalty for new blocks").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(XO_PENALTY_FLAG, Integer.class, CommonFlags.INT, "override penalty for crossovers").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    mNewPenalty = mFlags.isSet(NEW_PENALTY_FLAG) ? (int) mFlags.getValue(NEW_PENALTY_FLAG) : Search.DEFAULT_NEW_PENALTY;
    mXoPenalty = mFlags.isSet(XO_PENALTY_FLAG) ? (int) mFlags.getValue(XO_PENALTY_FLAG) : Search.DEFAULT_XO_PENALTY;
    final File genomeFile = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    try (final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(genomeFile)) {
      mPloidyMap = constructPloidyMap(sr);
    }
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final File bedFile = FileUtils.getZippedFileName(gzip, (File) mFlags.getValue(OUTPUT_REGIONS_FLAG));
    mSearch = new Search(mNewPenalty, mXoPenalty);
    try (final PrintStream outStream = new PrintStream(FileUtils.getZippedFileName(gzip, (File) mFlags.getValue(OUTPUT_FLAG)))) {
      mOut = outStream;
      try (final BedWriter bedOut = new BedWriter(FileUtils.createOutputStream(bedFile))) {
        mBed = bedOut;

        try (final VcfReader reader = VcfReader.openVcfReader((File) mFlags.getValue(VCF_FLAG))) {
          final VcfHeader header = reader.getHeader();
          for (final ContigField line : header.getContigLines()) {
            mChrLengths.put(line.getId(), line.getLength());
          }
          final Map<String, Sex> sampleSexMap = new HashMap<>();
          for (final SampleField line : header.getSampleLines()) {
            sampleSexMap.put(line.getId(), line.getSex());
          }
          final Map<String, Integer> sampleIndexMap = new HashMap<>();
          final List<String> sampleNames = header.getSampleNames();
          for (int i = 0; i < sampleNames.size(); ++i) {
            sampleIndexMap.put(sampleNames.get(i), i);
          }
          final String fatherName = (String) mFlags.getValue(FATHER_FLAG);
          final String motherName = (String) mFlags.getValue(MOTHER_FLAG);
          final int sampleFather = sampleIndexMap.get(fatherName);
          final int sampleMother = sampleIndexMap.get(motherName);
          final Sex[] sampleSex = getSampleSex(header.getSampleNames(), sampleSexMap, sampleFather, sampleMother);
          final Map<String, Integer> mismatchPloidyCounts = new LinkedHashMap<>();
          while (reader.hasNext()) {
            final VcfRecord rec = reader.next();
            final FamilyGt gty;
            try {
              gty = getFamilyGt(rec, sampleFather, sampleMother, sampleSex, mPloidyMap);
            } catch (final MismatchingPloidyException e) {
              Integer val = mismatchPloidyCounts.get(rec.getSequenceName());
              if (val == null) {
                val = 0;
              }
              mismatchPloidyCounts.put(rec.getSequenceName(), val + 1);
              continue;
            }
            final String seq = gty.seq();
            if (!seq.equals(mLastSeq)) {
              end();
              mOut.println("NS " + seq);
              mLastSeq = seq;
              mLastLength = gty.length();
              mBlock = null;
            }
            if (gty.allPloidyNone()) {
              mOut.println("NO " + gty);
              continue;
            }
            final boolean mendelian = gty.isMendelian();
            final boolean parentsSingleAllele = gty.parentsSingleAllele();
            final boolean allhet = gty.isAllHeterozygous();
            final String prefix = mendelian ? (parentsSingleAllele ? "HO" : (allhet ? "AH" : "OK")) : "ME";
            mOut.println(prefix + " " + gty);
            if (mendelian && !parentsSingleAllele && !allhet) {
              if (mBlock == null) {
                mBlock = new SegregationBlock(gty);
              } else {
                if (!mBlock.extend(gty)) {
                  mSearch.add(mBlock);
                  mBlock = new SegregationBlock(gty);
                }
              }
              //mOut.println("** " + mBlock.toString()); //debugging
            }
          }
          end();
          for (final Entry<String, Integer> seq : mismatchPloidyCounts.entrySet()) {
            Diagnostic.warning("There were " + seq.getValue() + " variants with genotypes that did not match the expected ploidy in chromosome " + seq.getKey());
          }
        }
      }
    }
    if (GzipAsynchOutputStream.BGZIP && gzip && !mFlags.isSet(CommonFlags.NO_INDEX)) {
      try {
        new TabixIndexer(bedFile).saveBedIndex();
      } catch (final UnindexableDataException e) {
        Diagnostic.warning(TabixIndexer.getTabixWarningMessage(bedFile, e));
      }
    }
    return 0;
  }

  static Map<Pair<Sex, String>, ReferenceSequence> constructPloidyMap(final SequencesReader sr) throws IOException {
    final Map<Pair<Sex, String>, ReferenceSequence> ploidyMap = new HashMap<>();
    if (ReferenceGenome.hasReferenceFile(sr)) {
      for (final Sex s : Sex.values()) {
        final ReferenceGenome ref = new ReferenceGenome(sr, s);
        for (long i = 0; i < sr.numberSequences(); ++i) {
          final String name = sr.name(i);
          ploidyMap.put(new Pair<>(s, name), ref.sequence(name));
        }
      }
    }
    return ploidyMap;
  }

  private void end() throws IOException {
    if (mBlock != null) {
      mSearch.add(mBlock);
    }
    final Iterable<SearchContainer> bestResult = mSearch.bestResult();
    if (bestResult == null) {
      return;
    }
    final RegionIterator reg = new RegionIterator(mOut, mBed, mLastSeq, mLastLength, mChrLengths, mPloidyMap);
    final Iterator<SearchContainer> it = bestResult.iterator();
    //TODO deal with X chromosome
    reg.iterate(it);
    mSearch = new Search(mNewPenalty, mXoPenalty);
  }

  private FamilyGt getFamilyGt(VcfRecord rec, int sampleFather, int sampleMother, Sex[] sexes, Map<Pair<Sex, String>, ReferenceSequence> ploidyMap) throws MismatchingPloidyException {
    final String[] sortedGts = new String[rec.getNumberOfSamples()];
    final List<String> gts = rec.getFormat(VcfUtils.FORMAT_GENOTYPE);
    sortedGts[0] = gts.get(sampleFather);
    sortedGts[1] = gts.get(sampleMother);
    int index = 2;
    for (int i = 0; i < gts.size(); ++i) {
      if (i == sampleFather || i == sampleMother) {
        continue;
      }
      sortedGts[index] = gts.get(i);
      ++index;
    }
    return FamilyGt.familyPloidy(rec.getSequenceName(), rec.getOneBasedStart(), sortedGts, sexes, ploidyMap);
  }

  protected static Sex[] getSampleSex(List<String> sampleNames, Map<String, Sex> sampleSexMap, int sampleFather, int sampleMother) {
    final Sex[] sampleSex = new Sex[sampleNames.size()];
    sampleSex[0] = sampleSexMap.get(sampleNames.get(sampleFather));
    sampleSex[1] = sampleSexMap.get(sampleNames.get(sampleMother));
    int index = 2;
    for (int i = 0; i < sampleNames.size(); ++i) {
      if (i == sampleFather || i == sampleMother) {
        continue;
      }
      sampleSex[index] = sampleSexMap.get(sampleNames.get(i));
      ++index;
    }
    return sampleSex;
  }

}
