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
package com.rtg.variant.cnv.segment;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INPUT_FLAG;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.bed.BedRangeLoader;
import com.rtg.bed.BedUtils;
import com.rtg.bed.NamedBedRangeLoader;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;

/**
 * Provide construction of a panel sample for CNV calling.
 */
public class CnvSummaryCli extends AbstractCli {

  private static final String SUMMARY_FLAG = "summary-regions";
  private static final String ALL_FLAG = "all-regions";

  @Override
  public String moduleName() {
    return "cnvsummary";
  }

  @Override
  public String description() {
    return "summarize the intersection of CNVs with regions of interest";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Summarize the intersection of CNVs with regions of interest.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('i', INPUT_FLAG, File.class, FILE, "VCF file containing CNV variants to be summarized. Use '-' to read from standard input").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, FILE, "BED output file. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired(SUMMARY_FLAG, File.class, FILE, "BED file supplying gene-scale regions to report CNV interactions with").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(ALL_FLAG, "if set, also report no alteration regions").setCategory(REPORTING);
    CommonFlags.initForce(mFlags);
    CommonFlags.initNoGzip(mFlags);

    mFlags.setValidator(flags -> CommonFlags.validateInputFile(flags, INPUT_FLAG)
      && CommonFlags.validateInputFile(flags, SUMMARY_FLAG)
      && CommonFlags.validateOutputFile(flags, FileUtils.getOutputFileName((File) flags.getValue(OUTPUT_FLAG), !flags.isSet(NO_GZIP), BedUtils.BED_SUFFIX))
    );
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final File vcfFile = (File) mFlags.getValue(INPUT_FLAG);
    final File outFile = (File) mFlags.getValue(OUTPUT_FLAG);
    final ReferenceRanges<String> reportRegions = BedRangeLoader.getReferenceRanges(new NamedBedRangeLoader(), (File) mFlags.getValue(SUMMARY_FLAG));
    final CnvSummaryReport reporter = new CnvSummaryReport(reportRegions, 0, !mFlags.isSet(ALL_FLAG));
    reporter.report(vcfFile, FileUtils.getOutputFileName(outFile, gzip, BedUtils.BED_SUFFIX));
    return 0;
  }
}
