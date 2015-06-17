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

import com.rtg.bed.BedUtils;
import com.rtg.calibrate.RecalibrateCli;
import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.ParamsTask;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamCommandHelper;
import com.rtg.usage.UsageMetric;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineOrientation;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * This is the top level module <code>map</code>, It contains all the user parameters with the correct name/syntax.
 */
public class MapCli extends ParamsCli<NgsParams>  {

  static final String MODULE_NAME = "map";
  static final String TOP_RANDOM = "top-random";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "read mapping";
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Aligns sequence reads onto a reference template, creating an alignments file in the Sequence Alignment/Map (SAM) format.");
    CommonFlagCategories.setCategories(flags);
    flags.setValidator(new MapFlagsValidator());
    MapFlags.initInputOutputFlags(flags);
    MapFlags.initPairedEndFormatFlags(flags);
    MapFlags.initPairedEndFlags(flags);
    MapFlags.initSharedFlags(flags);
    SamCommandHelper.initSamRg(flags);
    MapFlags.initMapFlags(flags);
    MapFlags.initSamOutputFlag(flags);
    MapFlags.initDontUnifyFlag(flags);
    initRamMapFlags(flags);
    CommonFlags.initIndexFlags(flags);
    MapFlags.initReadFreqFlag(flags, -1); //disable read frequency blocking
    MapFlags.initNoCalibrationFlag(flags);
    MapFlags.initSvPrepFlag(flags);
    MapFlags.initAlignerPenaltyFlags(flags);
    RecalibrateCli.bedFileFlag(flags);
  }

  @Override
  protected NgsParams makeParams() throws InvalidParamsException, IOException {
    return makeRamMapParams(mFlags, MapFlags.DEFAULT_WORD_SIZE, 1);
  }

  private static NgsParams makeRamMapParams(CFlags flags, int defWordSize, int defStepRatio) throws InvalidParamsException, IOException {
    final NgsParamsBuilder ngsParamsBuilder = NgsParams.builder();
    ngsParamsBuilder.useTopRandom(flags.isSet(TOP_RANDOM));
    final boolean paired = MapParamsHelper.isPaired(flags); //ReaderUtils.isPairedEndDirectory(input);
    final SAMReadGroupRecord rg;

    rg = MapParamsHelper.getSAMReadGroupRecord(flags);

    if (rg != null && (rg.getPlatform() == null || rg.getPlatform().length() == 0)) {
      Diagnostic.warning("Platform not specified in read group, it is recommended to set the platform tag.");
    }
    final NgsFilterParams filterParams = makeFilterParams(flags, paired, rg);
    final NgsOutputParams outputParams = makeOutputParams(flags, filterParams, rg);

    ngsParamsBuilder.outputParams(outputParams);
    if (flags.isSet(MapFlags.MIN_HITS_FLAG)) {
      ngsParamsBuilder.minHits((Integer) flags.getValue(MapFlags.MIN_HITS_FLAG));
    }

    ngsParamsBuilder.pairOrientation((MachineOrientation) flags.getValue(CommonFlags.PAIR_ORIENTATION_FLAG));
    ngsParamsBuilder.maxFragmentLength((Integer) flags.getValue(CommonFlags.MAX_FRAGMENT_SIZE));
    ngsParamsBuilder.minFragmentLength((Integer) flags.getValue(CommonFlags.MIN_FRAGMENT_SIZE));
    ngsParamsBuilder.compressHashes((Boolean) flags.getValue(MapFlags.COMPRESS_HASHES_FLAG));

    MapParamsHelper.populateAlignerPenaltiesParams(ngsParamsBuilder, flags);

    final long maxReadLength = MapParamsHelper.populateCommonMapParams(ngsParamsBuilder, flags, defWordSize, defStepRatio, flags.isSet(MapFlags.OUTPUT_READ_NAMES_FLAG), false);
    final NgsMaskParams maskParams = MapParamsHelper.makeMaskParams(flags, (int) maxReadLength, ngsParamsBuilder.mUseLongReadMapping, defWordSize);

    ngsParamsBuilder.maskParams(maskParams);
    return MapParamsHelper.createAndValidate(ngsParamsBuilder);
  }


  private static NgsOutputParams makeOutputParams(CFlags flags, NgsFilterParams filterParams, SAMReadGroupRecord rg) throws InvalidParamsException, IOException {
    final NgsOutputParamsBuilder ngsOutputParamsBuilder = NgsOutputParams.builder();
    ngsOutputParamsBuilder.progress(flags.isSet(BuildCommon.PROGRESS_FLAG))
    .outputDir((File) flags.getValue(CommonFlags.OUTPUT_FLAG))
    .tempFilesDir((File) flags.getValue(CommonFlags.TEMP_DIR))
    .filterParams(filterParams)
    .outputUnmated(!flags.isSet(MapFlags.NO_UNMATED))
    .outputUnmapped(!flags.isSet(MapFlags.NO_UNMAPPED))
    .tabular(false)
    .bam(!flags.isSet(MapFlags.SAM_FLAG))
    .unify(!flags.isSet(MapFlags.DONT_UNIFY_FLAG))
    .outputReadNames(flags.isSet(MapFlags.OUTPUT_READ_NAMES_FLAG))
    .outputIndex(!flags.isSet(CommonFlags.NO_INDEX));

    if (rg != null) {
      ngsOutputParamsBuilder.readGroup(rg);
    }

    ngsOutputParamsBuilder.calibrate(false);
    if (!flags.isSet(MapFlags.NO_CALIBRATION)) {
      if (flags.isSet(MapFlags.OUTPUT_UNFILTERED)) { //if all hits is set, do not calibrate
        Diagnostic.warning(MapFlags.OUTPUT_UNFILTERED + " is set, quality calibration output is disabled.");
      } else {
        if (rg == null) {
          Diagnostic.warning("No read group specified, quality calibration output is disabled.");
        } else {
          ngsOutputParamsBuilder.calibrate(true);
          if (flags.isSet(CommonFlags.BED_REGIONS_FLAG)) {
            ngsOutputParamsBuilder.calibrateRegions(BedUtils.regions((File) flags.getValue(CommonFlags.BED_REGIONS_FLAG)));
          }
        }
      }
    }
    ngsOutputParamsBuilder.svprep(false);
    if (MapParamsHelper.isPaired(flags) && !flags.isSet(MapFlags.NO_SVPREP)) {
      if (flags.isSet(MapFlags.OUTPUT_UNFILTERED)) { //if all hits is set, do not svprep
        Diagnostic.warning(MapFlags.OUTPUT_UNFILTERED + " is set, svprep output is disabled.");
      } else if (rg == null) {
        Diagnostic.warning("No read group specified, svprep output is disabled.");
      } else if (rg.getPlatform() == null || rg.getPlatform().length() == 0) {
        Diagnostic.warning("No platform in specified read group, svprep output is disabled.");
      } else if (ReadGroupUtils.platformToMachineType(rg, true) == null || ReadGroupUtils.platformToMachineType(rg, true).orientation() == null) {
        Diagnostic.warning("Unsupported platform in specified read group, svprep output is disabled.");
      } else {
        ngsOutputParamsBuilder.svprep(true);
      }
    }

    return ngsOutputParamsBuilder.create();
  }

  static NgsFilterParams makeFilterParams(CFlags flags, boolean paired, SAMReadGroupRecord rg) {
    final NgsFilterParams.NgsFilterParamsBuilder ngsFilterParamsBuilder = NgsFilterParams.builder();
    final OutputFilter filter;
    if (flags.isSet(MapFlags.OUTPUT_NULLFILTERED)) {
      filter = OutputFilter.NULL;
    } else if (flags.isSet(MapFlags.OUTPUT_UNFILTERED)) {
      filter = OutputFilter.SAM_UNFILTERED;
    } else if (paired) {
      filter = OutputFilter.TOPN_PAIRED_END;
    } else {
      filter = OutputFilter.SAM_SINGLE_END;
    }

    final int topN = (Integer) flags.getValue(MapFlags.TOPN_RESULTS_FLAG);
    final int maxTopResults = (Integer) flags.getValue(MapFlags.MAX_TOP_RESULTS_FLAG);
    ngsFilterParamsBuilder.outputFilter(filter)
    .zip(!flags.isSet(CommonFlags.NO_GZIP))
    .topN(Math.max(topN, maxTopResults))
    .maxTopResults(maxTopResults)
    .errorLimit((Integer) flags.getValue(MapFlags.XSCORE_INDEL)) //this used to be default for solexa / illumina short reads
    .exclude(flags.isSet(CommonFlags.EXCLUDE_FLAG));
    MapParamsHelper.populateAlignmentScoreSettings(flags, ngsFilterParamsBuilder, paired, rg);
    return ngsFilterParamsBuilder.create();
  }

  static void initRamMapFlags(CFlags flags) {
    flags.registerOptional('n', MapFlags.MAX_TOP_RESULTS_FLAG, Integer.class, CommonFlags.INT, "maximum number of top equal results output per read", 5).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.NO_UNMAPPED, "do not output unmapped reads").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.NO_UNMATED, "do not output unmated reads when in paired-end mode").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.OUTPUT_UNFILTERED, "output all alignments meeting thresholds instead of applying mating and N limits").setCategory(CommonFlagCategories.REPORTING);

    //x flags
    flags.registerOptional(MapFlags.X_LONG_READ, "use the non-default version for long read").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(MapFlags.MASK_FLAG, String.class, "STRING", "mask class name").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MapFlags.COMPRESS_HASHES_FLAG, Boolean.class, "BOOL", "compress hashes in indexes", true).setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(MapFlags.XSCORE_INDEL, Integer.class, CommonFlags.INT, "set max score indel for topn threshold", MapFlags.MAX_SCORE).setCategory(CommonFlagCategories.REPORTING); //7 was used for illumina mappings
    flags.registerOptional(MapFlags.OUTPUT_NULLFILTERED, "write nothing").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.MIN_HITS_FLAG, Integer.class, CommonFlags.INT, "Require this many hits to a logical read position before further processing").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(CommonFlags.EXCLUDE_FLAG, BuildCommon.RESOURCE.getString("EXCLUDE_DESC")).setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(MapFlags.TOPN_RESULTS_FLAG, Integer.class, CommonFlags.INT, "set the number of results per read for topn. Allowed values are between 1 and 255", 5).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(TOP_RANDOM, "output a single random top hit per read").setCategory(CommonFlagCategories.REPORTING);
  }

  static class MapFlagsValidator implements Validator  {
    @Override
    public boolean isValid(final CFlags flags) {
      return flags.checkNand(MapCli.TOP_RANDOM, MapFlags.OUTPUT_UNFILTERED) && validateInputOutput(flags) && validateParams(flags);
    }

    private boolean validateParams(final CFlags flags) {
      return MapFlags.validateMapParams(flags) && validateRamMapParams(flags);
    }

    private boolean validateInputOutput(final CFlags flags) {
      return MapFlags.validateMapInputOutputParams(flags) && validateRamMapInputOutputParams(flags);
    }

    private static boolean validateRamMapParams(CFlags flags) {
      if (!CommonFlags.validateFlagBetweenValues(flags, MapFlags.TOPN_RESULTS_FLAG, 1, 255)) {
        return false;
      }
      if (!CommonFlags.validateFlagBetweenValues(flags, MapFlags.MAX_TOP_RESULTS_FLAG, 1, 65535)) {
        return false;
      }
      if (flags.isSet(MapFlags.XSCORE_INDEL)) {
        if (!CommonFlags.validateFlagBetweenValues(flags, MapFlags.XSCORE_INDEL, 0, MapFlags.MAX_SCORE)) {
          return false;
        }
      }
      if (!MapFlags.validateSexTemplateReference(flags)) {
        return false;
      }
      if (!RecalibrateCli.checkBedFileFlag(flags)) {
        return false;
      }
      return true;
    }

    private static boolean validateRamMapInputOutputParams(CFlags flags) {
      return !(flags.isSet(SamCommandHelper.SAM_RG) && !SamCommandHelper.validateSamRg(flags));
    }

  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected ParamsTask<?, ?> task(final NgsParams params, final OutputStream out) {
    Diagnostic.userLog(MODULE_NAME + " paired=" + params.paired() + ", long=" + params.useLongReadMapping() + " running NgsTask");
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new NgsTask(params, out, usageMetric);
  }

}
