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

import static com.rtg.launcher.CommonFlags.TEMP_DIR;
import static com.rtg.ngs.MapFlags.COMPRESS_HASHES_FLAG;
import static com.rtg.ngs.MapFlags.TEMP_FILES_COMPRESSED;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.HashSet;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.bed.BedUtils;
import com.rtg.index.hash.ngs.FactoryUtil;
import com.rtg.index.hash.ngs.instances.AbstractCG2Mask;
import com.rtg.index.hash.ngs.instances.AbstractCGMask;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.DefaultReaderParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.ParamsTask;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.AlternatingSequencesWriter;
import com.rtg.reader.CgUtils;
import com.rtg.reader.CompressedMemorySequencesReader;
import com.rtg.reader.FormatCli;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SdfUtils;
import com.rtg.reader.SourceFormat;
import com.rtg.reader.TsvSequenceDataSource;
import com.rtg.reference.Sex;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamCommandHelper;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.License;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Top level module to run CG mapping
 */
public class CgMapCli extends ParamsCli<NgsParams> {

  private static final int MAX_INSERT_SIZE = 1000000000;

  private static final String MODULE_NAME = "cgmap";

  private static final String MAX_TOPN_RESULTS = "Xmax-topn-results";
  private static final String XINTSET_WINDOW = "Xintset-window";
  private static final String MASK_FLAG = "mask";
  private static final String OUTPUT_UNFILTERED = "all-hits";
  private static final String LEGACY_CIGARS = "legacy-cigars";

  private static final Validator VALIDATOR = new CgMapFlagsValidator();

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "read mapping for Complete Genomics data";
  }

  @TestClass("com.rtg.ngs.CgMapFlagsValidatorTest")
  private static class CgMapFlagsValidator implements Validator  {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      final boolean sdf = FormatCli.SDF_FORMAT.equals(flags.getValue(FormatCli.FORMAT_FLAG));
      if (!CommonFlags.validateReads(flags, sdf) || !CommonFlags.validateTemplate(flags)) {
        return false;
      }
      final File inputFile = (File) flags.getFlag(CommonFlags.READS_FLAG).getValue();
      if (sdf && !ReaderUtils.isPairedEndDirectory(inputFile)) {
        Diagnostic.error(ErrorType.NOT_A_PAIRED_END_SDF, inputFile.getPath());
        return false;
      }
      final String mask = flags.getValue(MASK_FLAG).toString();
      if (!FactoryUtil.checkMaskExists(mask)) {
        Diagnostic.error(ErrorType.INVALID_MASK, mask);
        return false;
      }

      return MapFlags.validateSexTemplateReference(flags)
      && flags.checkInRange(MapFlags.XSCORE_INDEL, 0, MapFlags.MAX_SCORE)
      && MapFlags.validateMinMaxFragmentSize(flags)
      && SamCommandHelper.validateSamRg(flags)
      && MapFlags.checkPercentRepeatFrequency(flags)
      && CommonFlags.validateThreads(flags)
      && flags.checkInRange(MapFlags.MAX_TOP_RESULTS_FLAG, 1, 255)
      && flags.checkNand(MapFlags.BAM_FLAG, CommonFlags.NO_GZIP)
      && CommonFlags.validateInputFile(flags, CommonFlags.BED_REGIONS_FLAG);
    }
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Aligns Complete Genomics sequence reads onto a reference template, creating an alignments file in the Sequence Alignment/Map (SAM) format.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setValidator(VALIDATOR);

    mFlags.registerRequired('i', CommonFlags.READS_FLAG, File.class, "SDF|FILE", "the Complete Genomics read set").setCategory(INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(mFlags);
    mFlags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF containing template to map against").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(MapFlags.NO_INMEMORY_TEMPLATE, "do not load the template in memory").setCategory(UTILITY);

    final Flag<String> maskFlag = mFlags.registerRequired(MASK_FLAG, String.class, CommonFlags.STRING, "read indexing method").setCategory(SENSITIVITY_TUNING);
    if (License.isDeveloper()) {
      maskFlag.setParameterRange(new String[]{
        "cg1", "cg1-fast", "cg2",
        "cgmaska15b1", "cgmaska1b1", "cgmaska15b1alt", "cgmaska1b1alt",
        "cg2maska1", "cg2maska11", "cg2maska15", "cg2maskw18", "cg2maskw18a1"
      });
    } else {
      maskFlag.setParameterRange(new String[]{"cg1", "cg1-fast", "cg2"});
    }
    MapFlags.initPairedEndFlags(mFlags);
    MapFlags.initSharedFlagsOnly(mFlags, IntegerOrPercentage.valueOf("95%"), 1, 1000);

    mFlags.registerOptional('e', MapFlags.MATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches allowed for mated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_MATED_MISMATCH_THRESHOLD)).setCategory(REPORTING);
    mFlags.registerOptional('E', MapFlags.UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches allowed for unmated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_UNMATED_MISMATCH_THRESHOLD)).setCategory(REPORTING);
    mFlags.registerOptional(MapFlags.NO_UNMATED, "do not report unmated").setCategory(UTILITY);
    mFlags.registerOptional(MapFlags.NO_UNMAPPED, "do not report unmapped").setCategory(UTILITY);

    mFlags.registerOptional(TEMP_DIR, File.class, CommonFlags.DIR, "directory used for temporary files (Defaults to output directory)").setCategory(UTILITY);

    mFlags.registerOptional('n', MapFlags.MAX_TOP_RESULTS_FLAG, Integer.class, CommonFlags.INT, "maximum number of top equal results output per read", 5).setCategory(REPORTING);
    final Flag<String> format = mFlags.registerOptional('F', FormatCli.FORMAT_FLAG, String.class, "FORMAT", "format of read data", FormatCli.SDF_FORMAT).setCategory(INPUT_OUTPUT);
    format.setParameterRange(new String[]{FormatCli.SDF_FORMAT, FormatCli.TSV_FORMAT});
    mFlags.registerOptional(LEGACY_CIGARS, "use legacy cigars in output").setCategory(UTILITY);

    CommonFlags.initReadRange(mFlags);

    SamCommandHelper.initSamRg(mFlags, "COMPLETE", REPORTING);
    CommonFlags.initIndexFlags(mFlags);
    //--X flags
    MapFlags.initReadFreqFlag(mFlags, 65535); //disable read frequency blocking
    mFlags.registerOptional(TEMP_FILES_COMPRESSED, Boolean.class, "BOOL", "gzip temporary SAM files", Boolean.TRUE).setCategory(REPORTING);
    mFlags.registerOptional(MapFlags.XSCORE_INDEL, Integer.class, CommonFlags.INT, "maximum score indel threshold", MapFlags.MAX_SCORE).setCategory(REPORTING);
    mFlags.registerOptional(MAX_TOPN_RESULTS, Integer.class, CommonFlags.INT, "sets the number of results per read for topn. Allowed values are between 1 and 255", 5).setCategory(REPORTING);
    mFlags.registerOptional(XINTSET_WINDOW, Integer.class, CommonFlags.INT, "windows for int set", 1).setCategory(UTILITY);
    mFlags.registerOptional(COMPRESS_HASHES_FLAG, Boolean.class, "BOOL", "compress hashes in indexes", Boolean.TRUE).setCategory(UTILITY);
    mFlags.registerOptional(OUTPUT_UNFILTERED, "output all alignments meeting thresholds instead of applying mating and N limits").setCategory(REPORTING);
    mFlags.registerOptional(MapFlags.N_AS_MISMATCH, "treat unknowns as mismatches").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MapFlags.SEX_FLAG, Sex.class, "sex", "sex of individual", null).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(CommonFlags.PEDIGREE_FLAG, File.class, CommonFlags.FILE, "genome relationships pedigree containing sex of sample").setCategory(SENSITIVITY_TUNING);
    MapFlags.initSamOutputFlag(mFlags);
    MapFlags.initDontUnifyFlag(mFlags);
    MapFlags.initNoCalibrationFlag(mFlags);
    MapFlags.initSvPrepFlag(mFlags);

  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected NgsParams makeParams() throws IOException {

    final NgsFilterParams filterParams = makeFilterParams();
    final NgsOutputParams outputParams = makeOutputParams(filterParams);


    final String maskName = (String) mFlags.getValue(MASK_FLAG);
    final NgsMaskParams maskParams = new NgsMaskParamsExplicit(maskName);
    final NgsParamsBuilder ngsParamBuilder = NgsParams.builder();
    ngsParamBuilder.name(mFlags.getName());
    ngsParamBuilder.compressHashes((Boolean) mFlags.getValue(COMPRESS_HASHES_FLAG));
    final Collection<ListenerType> listeners = new HashSet<>();
    listeners.add(ListenerType.CLI);
    ngsParamBuilder.listeners(listeners)
    .numberThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)))
    .readFreqThreshold((Integer) mFlags.getValue(MapFlags.READ_FREQUENCY_FLAG))
    .parallelUnmatedProcessing((Boolean) mFlags.getValue(MapFlags.PARALLEL_UNMATED_PROCESSING_FLAG));

    MapParamsHelper.populateProportionalRepeat(mFlags, ngsParamBuilder);

    final File reads = (File) mFlags.getValue(CommonFlags.READS_FLAG);
    final LongRange buildReaderRestriction = CommonFlags.getReaderRestriction(mFlags);
    try {
      if (FormatCli.getFormat(mFlags, false).getSourceFormat() == SourceFormat.SDF) {
        ngsParamBuilder.buildFirstParams(SequenceParams.builder().directory(ReaderUtils.getLeftEnd(reads)).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).readerRestriction(buildReaderRestriction).create());
        ngsParamBuilder.buildSecondParams(SequenceParams.builder().directory(ReaderUtils.getRightEnd(reads)).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).readerRestriction(buildReaderRestriction).create());

      } else {
        final TsvSequenceDataSource tsv = new TsvSequenceDataSource(reads, 5);
        final AlternatingSequencesWriter asw = new AlternatingSequencesWriter(tsv, null, PrereadType.CG, true);
        asw.setSdfId(new SdfId(0));
        final CompressedMemorySequencesReader[] readers = asw.processSequencesInMemoryPaired(reads, true, null, null, buildReaderRestriction);
        ngsParamBuilder.buildFirstParams(SequenceParams.builder().readerParam(new DefaultReaderParams(readers[0], buildReaderRestriction, SequenceMode.UNIDIRECTIONAL)).useMemReader(true).mode(SequenceMode.UNIDIRECTIONAL).readerRestriction(buildReaderRestriction).create()); // Reads
        ngsParamBuilder.buildSecondParams(SequenceParams.builder().readerParam(new DefaultReaderParams(readers[1], buildReaderRestriction, SequenceMode.UNIDIRECTIONAL)).useMemReader(true).mode(SequenceMode.UNIDIRECTIONAL).readerRestriction(buildReaderRestriction).create()); // Reads
      }
    } catch (final IOException e) {
      throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, reads.getPath());
    }
    final File template = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    final boolean inMemoryTemplate = !mFlags.isSet(MapFlags.NO_INMEMORY_TEMPLATE);
    SdfUtils.validateHasNames(template);
    ngsParamBuilder.outputParams(outputParams);
    final Sex sex = MapParamsHelper.getMappingSex(ngsParamBuilder, mFlags);
    final SequenceParams tParams = SequenceParams.builder().directory(template).mode(SequenceMode.UNIDIRECTIONAL).sex(sex).loadNames(true).useMemReader(inMemoryTemplate).create();
    if (outputParams.calibrateRegions() != null) {
      ReaderUtils.validateRegions(tParams.reader(), outputParams.calibrateRegions());
    }
    try {
      SdfUtils.validateNoDuplicates(tParams.reader(), true);
      ngsParamBuilder.searchParams(tParams)
      .useLongReadMapping(false)
      .maskParams(maskParams)
      .pairOrientation((MachineOrientation) mFlags.getValue(CommonFlags.PAIR_ORIENTATION_FLAG))
      .maxFragmentLength((Integer) mFlags.getValue(CommonFlags.MAX_FRAGMENT_SIZE))
      .minFragmentLength((Integer) mFlags.getValue(CommonFlags.MIN_FRAGMENT_SIZE))
      .intsetWindow((Integer) mFlags.getValue(XINTSET_WINDOW))
      .legacyCigars(mFlags.isSet(LEGACY_CIGARS));


      //TODO cg aligner sets some of these values internally - really should use these instead.
      // however the abstract temp file writer uses these values and so they need to be set here anyway.
      ngsParamBuilder.substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1);
      ngsParamBuilder.unknownsPenalty(mFlags.isSet(MapFlags.N_AS_MISMATCH) ? 1 : 0);  //cg only supports 1/1/1/0 or 1/1/1/1

      //cg sets up max shift as hard coded '7' in EditDistanceFactory, this is equivalent to that for 35 long reads
      ngsParamBuilder.alignerBandWidthFactor(new MaxShiftFactor(.5));

      final NgsParams params = ngsParamBuilder.create();
      final ISequenceParams leftParams = params.buildFirstParams();
      final ISequenceParams rightParams = params.buildSecondParams();
      if (leftParams.reader().getPrereadType() != PrereadType.CG) {
        throw new InvalidParamsException(ErrorType.NOT_A_CG_SDF, reads.toString());
      }
      if (rightParams.reader().getPrereadType() != PrereadType.CG) {
        throw new InvalidParamsException(ErrorType.NOT_A_CG_SDF, reads.toString());
      }
      if (leftParams.reader().numberSequences() != rightParams.reader().numberSequences()) {
        throw new InvalidParamsException("Left and right SDFs for read pair must have same number of sequences, actually had: "
            + leftParams.numberSequences() + " and " + rightParams.numberSequences());
      }
      final long length = leftParams.reader().maxLength();
      if (length != leftParams.reader().minLength()
        || length != rightParams.reader().minLength()
        || length != rightParams.reader().maxLength()) {
        throw new InvalidParamsException("Complete Genomics input must have all reads the same length.");
      }
      if (length != CgUtils.CG_RAW_READ_LENGTH
        && length != CgUtils.CG2_RAW_READ_LENGTH) {
        throw new InvalidParamsException("Complete Genomics input must have read length of " + CgUtils.CG_RAW_READ_LENGTH + " or " + CgUtils.CG2_RAW_READ_LENGTH + " bp");
      }
      if (length == CgUtils.CG_RAW_READ_LENGTH && !(maskParams.maskFactory(-1) instanceof AbstractCGMask.CGHashFunctionFactory)) {
        throw new InvalidParamsException("Mask '" + maskName + "' is not valid for CG version 1 reads, please select another via --mask");
      } else if (length == CgUtils.CG2_RAW_READ_LENGTH && !(maskParams.maskFactory(-1) instanceof AbstractCG2Mask.CG2HashFunctionFactory)) {
        throw new InvalidParamsException("Mask '" + maskName + "' is not valid for CG version 2 reads, please supply another via --mask");
      }
      final SAMReadGroupRecord samrg = params.outputParams().readGroup();
      if (samrg != null) {
        final String platform = samrg.getPlatform();
        final MachineType mt = length == CgUtils.CG_RAW_READ_LENGTH ? MachineType.COMPLETE_GENOMICS : MachineType.COMPLETE_GENOMICS_2;
        if (!mt.compatiblePlatform(platform)) {
          if (platform == null || platform.length() == 0) {
            Diagnostic.warning("Read group platform not set, defaulting to \"" + mt.platform() + "\"");
            samrg.setPlatform(mt.platform());
          } else {
            Diagnostic.warning("Read group platform is \"" + platform + "\", should be set to \"" + mt.platform() + "\"");
          }
        }
      }

      return params;
    } catch (final RuntimeException e) {
      tParams.close();
      throw e;
    }
  }

  private NgsOutputParams makeOutputParams(NgsFilterParams filterParams) throws IOException {
    final NgsOutputParamsBuilder ngsOutputParamsBuilder = NgsOutputParams.builder();
    ngsOutputParamsBuilder.outputIndex(!mFlags.isSet(CommonFlags.NO_INDEX))
    .progress(true)
    .outputDir((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG))
    .tempFilesDir((File) mFlags.getValue(TEMP_DIR))
    .filterParams(filterParams)
    .outputUnmated(!mFlags.isSet(MapFlags.NO_UNMATED))
    .outputUnmapped(!mFlags.isSet(MapFlags.NO_UNMAPPED))
    .sorted(false)
    .bam(!mFlags.isSet(MapFlags.SAM_FLAG))
    .unify(!mFlags.isSet(MapFlags.DONT_UNIFY_FLAG));

    final SAMReadGroupRecord rg = MapParamsHelper.getSAMReadGroupRecord(mFlags);
    if (rg != null) {
      ngsOutputParamsBuilder.readGroup(rg);
    }
    ngsOutputParamsBuilder.calibrate(false);
    if (!mFlags.isSet(MapFlags.NO_CALIBRATION)) {
      if (mFlags.isSet(OUTPUT_UNFILTERED)) { //if all hits is set, do not calibrate
        Diagnostic.warning(OUTPUT_UNFILTERED + " is set, quality calibration output is disabled.");
      } else {
        if (rg != null) {
          ngsOutputParamsBuilder.calibrate(true);
          if (mFlags.isSet(CommonFlags.BED_REGIONS_FLAG)) {
            ngsOutputParamsBuilder.calibrateRegions(BedUtils.regions((File) mFlags.getValue(CommonFlags.BED_REGIONS_FLAG)));
          }
        } else {
          Diagnostic.warning("No read group specified, quality calibration output is disabled.");
        }
      }
    }
    ngsOutputParamsBuilder.svprep(false);
    if (!mFlags.isSet(MapFlags.NO_SVPREP)) {
      if (mFlags.isSet(OUTPUT_UNFILTERED)) { //if all hits is set, do not svprep
        Diagnostic.warning(OUTPUT_UNFILTERED + " is set, svprep output is disabled.");
      } else if (rg == null) {
        Diagnostic.warning("No read group specified, svprep output is disabled.");
      } else if (ReadGroupUtils.platformToMachineType(rg, true) == null || ReadGroupUtils.platformToMachineType(rg, true).orientation() == null) {
        Diagnostic.warning("Unsupported platform in specified read group, svprep output is disabled.");
      } else {
        ngsOutputParamsBuilder.svprep(true);
      }
    }
    return ngsOutputParamsBuilder.create();
  }

  private NgsFilterParams makeFilterParams() {
    final NgsFilterParams.NgsFilterParamsBuilder ngsFilterParamsBuilder = NgsFilterParams.builder();
    final IntegerOrPercentage matedMismatches = (IntegerOrPercentage) mFlags.getValue(MapFlags.MATED_MISMATCH_THRESHOLD);
    final IntegerOrPercentage unmatedMismatches = (IntegerOrPercentage) mFlags.getValue(MapFlags.UNMATED_MISMATCH_THRESHOLD);
    final OutputFilter outputFilter;
    if (mFlags.isSet(OUTPUT_UNFILTERED)) {
      outputFilter = OutputFilter.SAM_UNFILTERED;
    } else {
      outputFilter = OutputFilter.TOPN_PAIRED_END;
    }
    final int topN = (Integer) mFlags.getValue(MapFlags.TOPN_RESULTS_FLAG);
    final int maxTopResults = (Integer) mFlags.getValue(MapFlags.MAX_TOP_RESULTS_FLAG);
    ngsFilterParamsBuilder.outputFilter(outputFilter)
    .zip(!mFlags.isSet(CommonFlags.NO_GZIP))
    .topN(Math.max(topN, maxTopResults))
    .maxTopResults(maxTopResults)
    .exclude(false)
    .useids(false)
    .matedMaxMismatches(matedMismatches)
    .unmatedMaxMismatches(unmatedMismatches)
    .errorLimit((Integer) mFlags.getValue(MapFlags.XSCORE_INDEL));

    return ngsFilterParamsBuilder.create();
  }

  @Override
  protected ParamsTask<?, ?> task(final NgsParams params, final OutputStream out) {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new NgsTask(params, out, usageMetric);
  }
}
