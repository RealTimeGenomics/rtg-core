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
import com.rtg.util.InvalidParamsException;
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
  protected NgsParams makeParams() throws InvalidParamsException, IOException {
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
    flags.registerExtendedHelp();
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

  private static NgsParams makeMapfParams(CFlags flags, int defWordSize, int defStepRatio) throws InvalidParamsException, IOException {
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
