/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
import com.rtg.util.intervals.RangeMeta;
import com.rtg.util.intervals.SimpleRangeMeta;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.VcfWriterFactory;
import com.rtg.vcf.header.VcfHeader;

/**
 * A class to take a BED file containing regions where we are confident of the phasing pattern,
 * and a VCF file containing variants and add annotations indicating if the variant
 * is consistent with the known phasing for the mother or father at that point (if there is
 * a consistent phasing known).
 */
public class SegregationCheckerCli extends AbstractCli {

  private static final String MODULE_NAME = "phasingeval";

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

    CommonFlags.initReferenceTemplate(mFlags, true);
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
    final File genomeFile = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    final Map<Pair<Sex, String>, ReferenceSequence> ploidyMap;
    try (final SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(genomeFile)) {
      ploidyMap = SegregationVcfSearch.constructPloidyMap(sr);
    }
    final File vcfOut = (File) mFlags.getValue(OUTPUT_FLAG);
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    try (final VcfReader reader = VcfReader.openVcfReader((File) mFlags.getValue(VCF_FLAG))) {
      final VcfHeader header = SegregationChecker.modifyHeader(reader.getHeader(), mFlags.isSet(REPAIR_FLAG));
      final File vcfFile = VcfUtils.getZippedVcfFileName(gzip, vcfOut);
      try (final VcfWriter writer = new VcfWriterFactory(mFlags).addRunInfo(true).make(header, vcfFile)) {
        new SegregationChecker((String) mFlags.getValue(FATHER_FLAG), (String) mFlags.getValue(MOTHER_FLAG), reader, writer, patterns, ploidyMap, mFlags.isSet(REPAIR_FLAG)).run();
      }
    }
    return 0;
  }

  private Map<String, RangeList<PatternHolder>> loadBed() throws IOException {
    final Map<String, List<RangeMeta<PatternHolder>>> ranges = new HashMap<>();
    try (final BedReader reader = BedReader.openBedReader(null, (File) mFlags.getValue(BED_FLAG), 2)) {
      while (reader.hasNext()) {
        final BedRecord rec = reader.next();
        if (!rec.getAnnotations()[0].matches("^[01?]+$") || !rec.getAnnotations()[1].matches("^[01?]+$")) {
          continue;
        }
        final List<RangeMeta<PatternHolder>> chrList = ranges.computeIfAbsent(rec.getSequenceName(), k -> new ArrayList<>());
        final int start = rec.getStart();
        int end = rec.getEnd();
        if (end == start) {
          ++end;
        }
        chrList.add(new SimpleRangeMeta<>(start, end, PatternHolder.fromPatternStrings(rec.getAnnotations()[0], rec.getAnnotations()[1], rec.getAnnotations()[2])));
      }
    }
    final Map<String, RangeList<PatternHolder>> patterns = new HashMap<>(ranges.size());
    for (final Map.Entry<String, List<RangeMeta<PatternHolder>>> chrEntry : ranges.entrySet()) {
      final RangeList<PatternHolder> rangeList = new RangeList<>(chrEntry.getValue());
      patterns.put(chrEntry.getKey(), rangeList);
    }
    return patterns;
  }

}
