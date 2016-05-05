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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.util.Pair;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.RangeList.RangeData;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.VcfHeader;

/**
 * A class to take a BED file containing regions where we are confident of the phasing pattern,
 * and a VCF file containing variants and add annotations indicating if the variant
 * is consistent with the known phasing for the mother or father at that point (if there is
 * a consistent phasing known).
 */
public class SegregationCheckerCli extends AbstractCli {

  private static final String MODULE_NAME = "phasingeval";

  private static final String TEMPLATE_FLAG = "template";
  private static final String BED_FLAG = "bed";
  private static final String VCF_FLAG = "vcf";
  private static final String FATHER_FLAG = "father";
  private static final String MOTHER_FLAG = "mother";
  private static final String OUTPUT_FLAG = "output";
  private static final String REPAIR_FLAG = "repair";


  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "evaluate a VCF with respect to phasing patterns";
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Annotate a VCF with phasing segregation information from a BED file.");

    mFlags.registerRequired('t', TEMPLATE_FLAG, File.class, "SDF", "SDF of the reference genome the reads have been mapped against").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(VCF_FLAG, File.class, CommonFlags.FILE, "input VCF file to be annotated").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(BED_FLAG, File.class, CommonFlags.FILE, "input BED file containing regions of phasing segregation information").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.FILE, "output VCF file containing new annotations").setCategory(INPUT_OUTPUT);

    mFlags.registerRequired(FATHER_FLAG, String.class, CommonFlags.STRING, "sample name of the father").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(MOTHER_FLAG, String.class, CommonFlags.STRING, "sample name of the mother").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(REPAIR_FLAG, "if set, repair variants where changing only one family member GT would allow consistency with phasing").setCategory(INPUT_OUTPUT);

    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final Map<String, RangeList<PatternHolder>> patterns = loadBed();
    final File genomeFile = (File) mFlags.getValue(TEMPLATE_FLAG);
    final Map<Pair<Sex, String>, ReferenceSequence> ploidyMap;
    try (final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(genomeFile)) {
      ploidyMap = SegregationVcfSearch.constructPloidyMap(sr);
    }
    final File vcfOut = (File) mFlags.getValue(OUTPUT_FLAG);
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);
    final boolean stdout = CommonFlags.isStdio(vcfOut);
    try (final VcfReader reader = VcfReader.openVcfReader((File) mFlags.getValue(VCF_FLAG))) {
      final VcfHeader header = SegregationChecker.modifyHeader(reader.getHeader(), mFlags.isSet(REPAIR_FLAG));
      final File vcfFile = stdout ? null : VcfUtils.getZippedVcfFileName(gzip, vcfOut);
      try (VcfWriter writer = new VcfWriter(header, vcfFile, out, gzip, index)) {
        new SegregationChecker((String) mFlags.getValue(FATHER_FLAG), (String) mFlags.getValue(MOTHER_FLAG), reader, writer, patterns, ploidyMap, mFlags.isSet(REPAIR_FLAG)).run();
      }
    }
    return 0;
  }

  private Map<String, RangeList<PatternHolder>> loadBed() throws IOException {
    final Map<String, List<RangeData<PatternHolder>>> ranges = new HashMap<>();
    try (final BedReader reader = BedReader.openBedReader(null, (File) mFlags.getValue(BED_FLAG), 2)) {
      while (reader.hasNext()) {
        final BedRecord rec = reader.next();
        if (!rec.getAnnotations()[0].matches("^[01?]+$") || !rec.getAnnotations()[1].matches("^[01?]+$")) {
          continue;
        }
        List<RangeData<PatternHolder>> chrList = ranges.get(rec.getSequenceName());
        if (chrList == null) {
          chrList = new ArrayList<>();
          ranges.put(rec.getSequenceName(), chrList);
        }
        final int start = rec.getStart();
        int end = rec.getEnd();
        if (end == start) {
          end++;
        }
        chrList.add(new RangeData<>(start, end, PatternHolder.fromPatternStrings(rec.getAnnotations()[0], rec.getAnnotations()[1], rec.getAnnotations()[2])));
      }
    }
    final Map<String, RangeList<PatternHolder>> patterns = new HashMap<>();
    for (final Map.Entry<String, List<RangeData<PatternHolder>>> chrEntry : ranges.entrySet()) {
      final RangeList<PatternHolder> rangeList = new RangeList<>(chrEntry.getValue());
      patterns.put(chrEntry.getKey(), rangeList);
    }
    return patterns;
  }

}
