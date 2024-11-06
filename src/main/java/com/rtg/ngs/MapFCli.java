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
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderUtils;
import com.rtg.sam.SamCommandHelper;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.machine.MachineOrientation;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Map filter CLI class
 */
public class MapFCli extends ParamsCli<NgsParams> {

  private static final int DEFAULT_WORD_SIZE = 22;

  @Override
  public String moduleName() {
    return "mapf";
  }

  @Override
  public String description() {
    return "read mapping for filtering purposes";
  }

  protected static class MapFTask extends NgsTask {

    public MapFTask(NgsParams params, OutputStream defaultOutput, UsageMetric metric) {
      super(params, defaultOutput, metric);
      if (params.paired()) {
        mStatistics = new MapFilterPairedMapStatistics(params.directory());
      } else {
        mStatistics = new MapFilterSingleEndMapStatistics(params.directory());
      }
    }
  }

  @Override
  protected NgsParams makeParams() throws IOException {
    return makeMapfParams(mFlags, DEFAULT_WORD_SIZE, 2);
  }

  @Override
  protected IORunnable task(NgsParams params, OutputStream out) {
    return new MapFTask(params, out, mUsageMetric);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  protected static void initFlags(CFlags flags) {
    flags.setDescription("Filters reads for contaminant sequences by mapping them against the contaminant reference. It outputs two SDF files, one containing the input reads that map to the reference and one containing those that do not.");
    CommonFlagCategories.setCategories(flags);
    flags.setValidator(new MapfFlagsValidator());
    //CommonFlags.initMapIOFlags(flags);
    MapFlags.initInputOutputFlags(flags);
    MapFlags.initPairedEndFormatFlags(flags);
    MapFlags.initSharedFlagsOnly(flags);
    MapFlags.initMaskFlagsOnly(flags);
    MapFlags.initWordSize(flags, "word size (Default is " + DEFAULT_WORD_SIZE + ")");
    MapFlags.initStepSize(flags, "step size (Default is half word size)");
    MapFlags.initReadFreqFlag(flags, 1);
    MapFlags.initMapFlags(flags);
    MapFlags.initBamOutputFlag(flags, "output the alignment files in BAM format");
    MapFlags.initSamOutputFlag(flags);
    MapFlags.initDontUnifyFlag(flags);
    CommonFlags.initIndexFlags(flags);
    MapFlags.initPairedEndFlags(flags);
    MapFlags.initAlignerPenaltyFlags(flags);
    SamCommandHelper.initSamRg(flags);
  }

  static class MapfFlagsValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (flags.isSet(MapFlags.DONT_UNIFY_FLAG)) {
        if (!(flags.isSet(MapFlags.SAM_FLAG) || flags.isSet(MapFlags.BAM_FLAG))) {
          flags.setParseMessage("Cannot use --" + MapFlags.DONT_UNIFY_FLAG + " unless --" + MapFlags.SAM_FLAG + " or --" + MapFlags.BAM_FLAG + " is set");
          return false;
        }
      }
      return validateInputOutput(flags) && validateParams(flags) && SamCommandHelper.validateSamRg(flags);
    }
    private boolean validateParams(final CFlags flags) {
      return MapFlags.validateMapParams(flags);
    }
    private boolean validateInputOutput(final CFlags flags) {
      return MapFlags.validateMapInputOutputParams(flags);
    }
  }

  private static NgsParams makeMapfParams(CFlags flags, int defWordSize, int defStepRatio) throws IOException {
    final NgsParamsBuilder ngsParamsBuilder = NgsParams.builder();

    final NgsFilterParams.NgsFilterParamsBuilder ngsFilterParamsBuilder = NgsFilterParams.builder();
    final OutputFilter filter = OutputFilter.SAM_UNFILTERED;
    final File input = (File) flags.getValue(CommonFlags.READS_FLAG);

    final SAMReadGroupRecord rg = MapParamsHelper.getSAMReadGroupRecord(flags);
    MapParamsHelper.populateAlignmentScoreSettings(flags, ngsFilterParamsBuilder, ReaderUtils.isPairedEndDirectory(input), rg);
    ngsFilterParamsBuilder.outputFilter(filter).zip(!flags.isSet(CommonFlags.NO_GZIP)).errorLimit(MapFlags.MAX_SCORE);
    final NgsFilterParams filterParams = ngsFilterParamsBuilder.create();

    final NgsOutputParams outputParams = makeOutputParamsMapf(flags, filterParams, rg);
    ngsParamsBuilder.outputParams(outputParams);

    final long readLength = MapParamsHelper.populateCommonMapParams(ngsParamsBuilder, flags, defWordSize, defStepRatio, true, true);

    final NgsMaskParams maskParams = MapParamsHelper.makeMaskParams(flags, (int) readLength, ngsParamsBuilder.mUseLongReadMapping, defWordSize);
    ngsParamsBuilder.maskParams(maskParams);
    ngsParamsBuilder.compressHashes(true);
    ngsParamsBuilder.pairOrientation((MachineOrientation) flags.getValue(CommonFlags.PAIR_ORIENTATION_FLAG));
    ngsParamsBuilder.maxFragmentLength((Integer) flags.getValue(CommonFlags.MAX_FRAGMENT_SIZE));
    ngsParamsBuilder.minFragmentLength((Integer) flags.getValue(CommonFlags.MIN_FRAGMENT_SIZE));

    MapParamsHelper.populateAlignerPenaltiesParams(ngsParamsBuilder, flags);

    return MapParamsHelper.createAndValidate(ngsParamsBuilder);
  }

  private static NgsOutputParams makeOutputParamsMapf(final CFlags flags, final NgsFilterParams filterParams, final SAMReadGroupRecord rg) {
    final NgsOutputParamsBuilder ngsOutputParamsBuilder = NgsOutputParams.builder();
    ngsOutputParamsBuilder.progress(flags.isSet(BuildCommon.PROGRESS_FLAG))
    .outputDir((File) flags.getValue(CommonFlags.OUTPUT_FLAG))
    .tempFilesDir((File) flags.getValue(CommonFlags.TEMP_DIR))
    .filterParams(filterParams)
    .outputUnmated(true)
    .outputUnmapped(true)
    .tabular(false)
    .outputIndex(!flags.isSet(CommonFlags.NO_INDEX))
    .bam(flags.isSet(MapFlags.BAM_FLAG))
    .unify(!flags.isSet(MapFlags.DONT_UNIFY_FLAG))
    .sam(flags.isSet(MapFlags.SAM_FLAG))
    .calibrate(false)
    .svprep(false)
    .sdf(true)
    .outputReadNames(flags.isSet(MapFlags.OUTPUT_READ_NAMES_FLAG));

    if (rg != null) {
      ngsOutputParamsBuilder.readGroup(rg);
    }

    return ngsOutputParamsBuilder.create();
  }

}
