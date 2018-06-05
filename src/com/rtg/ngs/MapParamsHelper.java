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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.TimeUnit;

import com.rtg.alignment.AlignerMode;
import com.rtg.index.BlacklistFilterMethod;
import com.rtg.index.FixedRepeatFrequencyFilterMethod;
import com.rtg.index.HashBlacklist;
import com.rtg.index.IndexFilterMethod;
import com.rtg.index.ProportionalRepeatFrequencyFilterMethod;
import com.rtg.index.UnionRepeatFrequencyFilterMethod;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.DefaultReaderParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.AlternatingSequencesWriter;
import com.rtg.reader.DataSourceDescription;
import com.rtg.reader.FormatCli;
import com.rtg.reader.IndexFile;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.RightSimpleNames;
import com.rtg.reader.SdfId;
import com.rtg.reader.SdfUtils;
import com.rtg.reader.SequenceDataSource;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesWriter;
import com.rtg.reader.SimpleNames;
import com.rtg.reader.SourceFormat;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SamCommandHelper;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.InputFileUtils;
import com.rtg.util.machine.MachineType;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Populate parameter objects for Map commands
 */
public final class MapParamsHelper {

  private MapParamsHelper() { }

  static long populateCommonMapParams(NgsParamsBuilder ngsParamsBuilder, CFlags flags, int defWordSize, int defStepRatio, boolean includeNames, boolean includeFullNames) throws IOException {
    final int numberThreads = CommonFlags.parseThreads((Integer) flags.getValue(CommonFlags.THREADS_FLAG));

    final NameParams nameParams = new NameParams(includeNames, includeFullNames);
    final long maxReadLength = populateReaders(ngsParamsBuilder, flags, nameParams);
    if (ngsParamsBuilder.mUseLongReadMapping) {
      MapParamsHelper.addLongReadParameters(flags, ngsParamsBuilder, maxReadLength, defWordSize, defStepRatio);
    }
    final int wordSize = MapFlags.getWordSize(flags, (int) maxReadLength, defWordSize);

    if (ngsParamsBuilder.mStepSize > maxReadLength) {
      throw new InvalidParamsException("Step size (" + ngsParamsBuilder.mStepSize + ") must be less than or equal to max read length (" + maxReadLength + ")");
    }

    final Flag<?> alignerChain = flags.getFlag(MapFlags.ALIGNER_MODE_FLAG);
    if (alignerChain != null) {
      ngsParamsBuilder.alignerMode((AlignerMode) alignerChain.getValue());
    }
    final Flag<?> singleIndelPenalties = flags.getFlag(MapFlags.SINGLE_INDEL_PENALTIES_FLAG);
    if (singleIndelPenalties != null) {
      ngsParamsBuilder.singleIndelPenalties((String) singleIndelPenalties.getValue());
    }

    final Collection<ListenerType> listeners = new HashSet<>();
    listeners.add(ListenerType.CLI);

    ngsParamsBuilder.listeners(listeners)
                    .numberThreads(numberThreads)
                    .threadMultiplier((Integer) flags.getValue(MapFlags.THREAD_MULTIPLIER))
                    .readFreqThreshold((Integer) flags.getValue(MapFlags.READ_FREQUENCY_FLAG))
                    .legacyCigars(flags.isSet(MapFlags.LEGACY_CIGARS))
                    .parallelUnmatedProcessing((Boolean) flags.getValue(MapFlags.PARALLEL_UNMATED_PROCESSING_FLAG));

    populateHashFilteringSettings(flags, ngsParamsBuilder, wordSize);
    return maxReadLength;
  }

  // Sets up appropriate repeat frequency filtering.
  // Defaults to read based repeat frequency filtering, but allows either (or both) to be explicitly requested
  static void populateHashFilteringSettings(CFlags flags, NgsParamsBuilder builder, int wordSize) throws IOException {
    final boolean forceBlacklist = flags.isSet(MapFlags.BLACKLIST_THRESHOLD);

    final boolean forceRepeat = flags.isSet(CommonFlags.REPEAT_FREQUENCY_FLAG)
      || flags.isSet(MapFlags.MAX_REPEAT_FREQUENCY_FLAG)
      || flags.isSet(MapFlags.MIN_REPEAT_FREQUENCY_FLAG);

    final List<IndexFilterMethod> filters = new ArrayList<>();
    if (forceBlacklist) {
      filters.add(makeReferenceRepeatFrequencyFilter(flags, builder, wordSize));
    }

    if (forceRepeat || !forceBlacklist) {
      filters.add(makeReadRepeatFrequencyFilter(flags));
    }
    assert filters.size() > 0;
    if (filters.size() == 1) {
      builder.indexFilter(filters.get(0));
    } else {
      builder.indexFilter(new UnionRepeatFrequencyFilterMethod(filters));
    }
    Diagnostic.userLog("Using IndexFilterMethod: " + builder.mIndexFilter);
  }

  /**
   * Populates the builder with repeat frequency related values. Note use <code>populateHashFilteringSettings</code> instead
   * if using the blacklist is an option
   * @param flags the flags
   * @param builder the builder
   */
  public static void populateReadRepeatFrequencyFilter(CFlags flags, NgsParamsBuilder builder) {
    builder.indexFilter(makeReadRepeatFrequencyFilter(flags));
  }

  private static IndexFilterMethod makeReferenceRepeatFrequencyFilter(CFlags flags, NgsParamsBuilder builder, int wordSize) throws IOException {
    if (!builder.mUseLongReadMapping) {
      throw new InvalidParamsException(ErrorType.BLACKLIST_LONG_READ_ONLY);
    }
    final File refDir = builder.searchParams().directory();
    final boolean blackListExists = HashBlacklist.blacklistExists(refDir, wordSize);
    if (!blackListExists) {
      throw new InvalidParamsException("A blacklist does not exist in " + refDir + " for word size " + wordSize);
    }
    final Integer blacklistThreshold = (Integer) flags.getValue(MapFlags.BLACKLIST_THRESHOLD);
    return BlacklistFilterMethod.loadBlacklist(refDir, wordSize, blacklistThreshold, builder.mNumberThreads);
  }

  private static IndexFilterMethod makeReadRepeatFrequencyFilter(CFlags flags) {
    final IndexFilterMethod filter;
    final IntegerOrPercentage repeat = (IntegerOrPercentage) flags.getValue(CommonFlags.REPEAT_FREQUENCY_FLAG);
    if (repeat.isPercentage()) {
      final int maxHC = (flags.getFlag(MapFlags.MAX_REPEAT_FREQUENCY_FLAG) != null) ? (Integer) flags.getValue(MapFlags.MAX_REPEAT_FREQUENCY_FLAG) : 1000;
      final int minHC = (flags.getFlag(MapFlags.MIN_REPEAT_FREQUENCY_FLAG) != null) ? (Integer) flags.getValue(MapFlags.MIN_REPEAT_FREQUENCY_FLAG) : 1;
      filter = new ProportionalRepeatFrequencyFilterMethod(100 - repeat.getValue(100), maxHC, minHC);
    } else {
      filter = new FixedRepeatFrequencyFilterMethod(repeat.getRawValue());
    }
    return filter;
  }

  static void populateAlignmentScoreSettings(CFlags flags, NgsFilterParams.NgsFilterParamsBuilder ngsFilterParamsBuilder, boolean paired, SAMReadGroupRecord rg)  {
    final boolean isIonTorrent = rg != null && MachineType.IONTORRENT.compatiblePlatform(rg.getPlatform());
    if (isIonTorrent) {
      Diagnostic.developerLog("IonTorrent mode enabled");
    }
    final IntegerOrPercentage ionTorrentDefaultThreshold = new IntegerOrPercentage("10%");
    if (paired) {
      final IntegerOrPercentage matedAS;
      if (flags.isSet(MapFlags.MATED_MISMATCH_THRESHOLD)) {
        matedAS = (IntegerOrPercentage) flags.getValue(MapFlags.MATED_MISMATCH_THRESHOLD);
      } else {
        if (!flags.isSet(MapFlags.MAX_ALIGNMENT_MISMATCHES) && isIonTorrent) {
          matedAS = ionTorrentDefaultThreshold;
        } else {
          matedAS = (IntegerOrPercentage) flags.getValue(MapFlags.MAX_ALIGNMENT_MISMATCHES);  //the defaults for mated alignment threshold and max align score are supposed to be the same
        }
      }
      final IntegerOrPercentage unmatedAS;
      if (flags.isSet(MapFlags.UNMATED_MISMATCH_THRESHOLD) || !isIonTorrent) { //if flag is set OR not iontorrent
        unmatedAS = (IntegerOrPercentage) flags.getValue(MapFlags.UNMATED_MISMATCH_THRESHOLD);
      } else {  // ie flag is not set AND this is iontorrent
        unmatedAS = ionTorrentDefaultThreshold;
      }
      if (flags.isSet(MapFlags.UNMATED_MISMATCH_THRESHOLD) && unmatedAS.compareTo(matedAS) > 0) {
        Diagnostic.warning("--" + MapFlags.UNMATED_MISMATCH_THRESHOLD + " should not be greater than --" + MapFlags.MATED_MISMATCH_THRESHOLD);
      }
      ngsFilterParamsBuilder.unmatedMaxMismatches(unmatedAS);
      ngsFilterParamsBuilder.matedMaxMismatches(matedAS);
    } else {
      if (!flags.isSet(MapFlags.MAX_ALIGNMENT_MISMATCHES) && isIonTorrent) {
        ngsFilterParamsBuilder.matedMaxMismatches(ionTorrentDefaultThreshold).unmatedMaxMismatches(ionTorrentDefaultThreshold);
      } else {
        ngsFilterParamsBuilder.matedMaxMismatches((IntegerOrPercentage) flags.getValue(MapFlags.MAX_ALIGNMENT_MISMATCHES))
              .unmatedMaxMismatches((IntegerOrPercentage) flags.getValue(MapFlags.MAX_ALIGNMENT_MISMATCHES));
      }
    }
  }

  static NgsParams createAndValidate(NgsParamsBuilder ngsParamsBuilder) throws IOException {
    final NgsParams localParams = ngsParamsBuilder.create();
    if (localParams.paired()) {
      if (localParams.buildFirstParams().reader().hasNames() != localParams.buildSecondParams().reader().hasNames()) {
        if (localParams.outputParams().outputReadNames()) {
          throw new InvalidParamsException("Names only present on one arms SDF and read names requested");
        }
      }
    }
    if (localParams.outputParams().outputReadNames() && !localParams.buildFirstParams().reader().hasNames()) {
      throw new InvalidParamsException("Names not present in SDF and read names requested");
    }
    if (localParams.outputParams().calibrateRegions() != null) {
      ReaderUtils.validateRegions(localParams.searchParams().reader(), localParams.outputParams().calibrateRegions());
    }

    localParams.globalIntegrity();
    return localParams;
  }

  static boolean isPaired(CFlags flags) {
    final DataSourceDescription inputDesc = FormatCli.getFormat(flags, false);
    final boolean paired;
    if (inputDesc.getSourceFormat() == SourceFormat.SDF) {
      final File reads = (File) flags.getValue(CommonFlags.READS_FLAG);
      paired = ReaderUtils.isPairedEndDirectory(reads);
    } else {
      paired = !flags.isSet(CommonFlags.READS_FLAG) || inputDesc.isInterleaved();
    }
    return paired;
  }

  /**
   *
   * @param ngsParamsBuilder builder to populate
   * @param flags flags to get configuration from
   * @return read length value
   * @throws InvalidParamsException when stuff goes wrong
   * @throws IOException when other stuff goes wrong
   */
  private static long populateReaders(NgsParamsBuilder ngsParamsBuilder, CFlags flags, NameParams nameParams) throws IOException {
    final boolean paired = initReaders(ngsParamsBuilder, flags, nameParams, true, SequenceMode.BIDIRECTIONAL, SequenceMode.UNIDIRECTIONAL, new SamSequenceReaderParams(false, true));
    final long maxReadLength;
    final boolean useLongReads;
    if (!paired) {
      final ISequenceParams params = ngsParamsBuilder.buildFirstParams();
      maxReadLength = params.reader().maxLength();
      final long minReadLength = params.reader().minLength();
      useLongReads = flags.isSet(MapFlags.FORCE_LONG_FLAG)
              || maxReadLength > 63
              || (minReadLength != maxReadLength);
      Diagnostic.userLog("Entering single end read mode read length=" + maxReadLength);
    } else {
      final ISequenceParams leftParams = ngsParamsBuilder.buildFirstParams();
      final ISequenceParams rightParams = ngsParamsBuilder.buildSecondParams();
      final long leftReadLength = leftParams.reader().maxLength();
      final long rightReadLength = rightParams.reader().maxLength();
      maxReadLength = Math.max(leftReadLength, rightReadLength);
      final long minReadLength = Math.min(leftParams.reader().minLength(), rightParams.reader().minLength());
      useLongReads = flags.isSet(MapFlags.FORCE_LONG_FLAG)
              || (maxReadLength > 63)
              || (minReadLength != maxReadLength);
      Diagnostic.userLog("Entering paired end read mode 1st arm read length=" + leftReadLength + " 2nd arm read length=" + rightReadLength);
    }

    if (maxReadLength == 0) {
      throw new InvalidParamsException(ErrorType.INFO_ERROR, "Read length must be greater than 0");
    }

    ngsParamsBuilder.useLongReadMapping(useLongReads);
    return maxReadLength;
  }

  /**
   * Initialise sequences readers
   * @param ngsParamsBuilder builder to populate
   * @param flags flags to get configuration from
   * @param nameParams options for read name loading
   * @param useQuality true to read quality when using FASTQ
   * @param templateMode the template sequence mode
   * @param readsMode the reads sequence mode
   * @param samParams parameters for when reading sequence data from SAM or BAM files
   * @return true if paired, false otherwise
   * @throws InvalidParamsException when stuff goes wrong
   * @throws IOException when other stuff goes wrong
   */
  public static boolean initReaders(NgsParamsBuilder ngsParamsBuilder, CFlags flags, NameParams nameParams, boolean useQuality, SequenceMode templateMode, SequenceMode readsMode, SamSequenceReaderParams samParams) throws IOException {
    final DataSourceDescription inputDesc = FormatCli.getFormat(flags, useQuality);
    final File reads = (File) flags.getValue(CommonFlags.READS_FLAG);
    final boolean paired = MapParamsHelper.isPaired(flags);
    if (!paired) {
      if (inputDesc.getSourceFormat() == SourceFormat.SDF) {
        makeSequenceParamsMulti(ngsParamsBuilder, flags, reads, null, nameParams, templateMode, readsMode);
      } else {
        makeOnTheFlySequenceParamsMulti(ngsParamsBuilder, flags, inputDesc, reads, null, nameParams, useQuality, templateMode, readsMode, samParams);
      }
      final ISequenceParams params = ngsParamsBuilder.buildFirstParams();
      if (params.reader().getPrereadType() == PrereadType.CG) {
        throw new InvalidParamsException(ErrorType.IS_A_CG_SDF, reads.getPath());
      }
    } else {
      if (inputDesc.getSourceFormat() == SourceFormat.SDF) {
        final File left = ReaderUtils.getLeftEnd(reads);
        final File right = ReaderUtils.getRightEnd(reads);
        makeSequenceParamsMulti(ngsParamsBuilder, flags, left, right, nameParams, templateMode, readsMode);
      } else if (inputDesc.isInterleaved()) {
        makeOnTheFlySequenceParamsMulti(ngsParamsBuilder, flags, inputDesc, reads, null, nameParams, useQuality, templateMode, readsMode, samParams);
      } else {
        final File leftFile = (File) flags.getValue(FormatCli.LEFT_FILE_FLAG);
        final File rightFile = (File) flags.getValue(FormatCli.RIGHT_FILE_FLAG);
        if (InputFileUtils.checkIdenticalPaths(leftFile, rightFile)) {
          throw new InvalidParamsException("Paths given for --" + FormatCli.LEFT_FILE_FLAG + " and --" + FormatCli.RIGHT_FILE_FLAG + " are the same file.");
        }
        makeOnTheFlySequenceParamsMulti(ngsParamsBuilder, flags, inputDesc, leftFile, rightFile, nameParams, useQuality, templateMode, readsMode, samParams);
      }
      final ISequenceParams leftParams = ngsParamsBuilder.buildFirstParams();
      final ISequenceParams rightParams = ngsParamsBuilder.buildSecondParams();
      if (leftParams.numberSequences() != rightParams.numberSequences()) {
        throw new InvalidParamsException("Left and right SDFs for read pair must have same number of sequences, actually had: "
                + leftParams.numberSequences() + " and " + rightParams.numberSequences());
      }
      if ((leftParams.reader().getPrereadType() == PrereadType.CG)
              || (rightParams.reader().getPrereadType() == PrereadType.CG)) {
        throw new InvalidParamsException(ErrorType.IS_A_CG_SDF, reads.getPath());
      }
    }
    return paired;
  }

  static NgsMaskParams makeMaskParams(final CFlags flags, final int readLength, boolean useLongReads, int defWordSize) {
    final NgsMaskParams maskParams;
    if (flags.isSet(MapFlags.MASK_FLAG)) {
      maskParams = new NgsMaskParamsExplicit((String) flags.getValue(MapFlags.MASK_FLAG));
    } else {
      final int w = MapFlags.getWordSize(flags, readLength, defWordSize);
      final int a = (Integer) flags.getValue(MapFlags.SUBSTITUTIONS_FLAG);
      final int b = (Integer) flags.getValue(CommonFlags.INDELS_FLAG);
      final int c = (Integer) flags.getValue(MapFlags.INDEL_LENGTH_FLAG);
      final int s = Math.max(a, b);
      if (readLength < w) {
        throw new InvalidParamsException(ErrorType.WORD_NOT_LESS_READ, w + "", readLength + "");
      }
      if (b > 0 && c > readLength - w) {
        throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-c", c + "", (readLength - w) + "");
      }
      if (b > 0 && c <= 0) {
        throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-c", c + "", "1");
      }
      if (b > readLength - w) {
        throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-b", b + "", (readLength - w) + "");
      }
      if (s > readLength - w) {
        throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-a", s + "", (readLength - w) + "");
      }
      if (useLongReads) {
        maskParams = new LongReadMaskParams(w, s, b, c);
      } else {
        maskParams = new NgsMaskParamsGeneral(w, s, b, c  /* no CG data here */);
      }
      if (!maskParams.isValid(readLength)) {
        throw new InvalidParamsException(ErrorType.INVALID_MASK_PARAMS);
      }
    }
    return maskParams;
  }

  private static void addLongReadParameters(final CFlags flags, final NgsParamsBuilder ngsParamsBuilder, final long readLength, final int defWordSize, int defStepRatio) {
    final int w = MapFlags.getWordSize(flags, (int) readLength, defWordSize);
    final int subs = (Integer) flags.getValue(MapFlags.SUBSTITUTIONS_FLAG);
    if (!((readLength / w) >= subs + 1)) {
      throw new InvalidParamsException(ErrorType.INVALID_LONG_READ_PARAMS, Integer.toString(w), Integer.toString(subs));
    }
    final int step;
    if (flags.isSet(MapFlags.STEP_FLAG)) {
      step = (Integer) flags.getValue(MapFlags.STEP_FLAG);
    } else {
      step = Math.max(w / defStepRatio, 1);
    }
    ngsParamsBuilder.stepSize(step);
  }


  private static FutureTask<SequenceParams> getTemplateFutureTask(NgsParamsBuilder ngsParamsBuilder, CFlags flags, boolean includeFullNames, SequenceMode templateMode) throws IOException {
    final boolean templateMem = !flags.isSet(MapFlags.NO_INMEMORY_TEMPLATE);
    final File template = (File) flags.getValue(CommonFlags.TEMPLATE_FLAG);
    final Sex sex = getMappingSex(ngsParamsBuilder, flags);
    return new FutureTask<>(new SequenceParamsCallableSdf(template, LongRange.NONE, templateMem, sex, includeFullNames, templateMode));
  }

  static Sex getMappingSex(NgsParamsBuilder ngsParamsBuilder, CFlags flags) throws IOException {
    final Sex sex;
    if (flags.isSet(CommonFlags.PEDIGREE_FLAG)) {
      final SAMReadGroupRecord rg = ngsParamsBuilder.mOutputParams.readGroup();
      if (rg == null) {
        throw new InvalidParamsException("No read group information has been provided, so cannot obtain sex information from pedigree file.");
      }
      final String sample = rg.getSample();
      if (sample == null) {
        throw new InvalidParamsException("Supplied read group information does not contain a sample, so cannot obtain sex information from pedigree file.");
      }
      final File relfile = (File) flags.getValue(CommonFlags.PEDIGREE_FLAG);
      final GenomeRelationships gr = GenomeRelationships.loadGenomeRelationships(relfile);
      if (!gr.hasGenome(sample)) {
        throw new InvalidParamsException("Supplied pedigree file does not contain sample " + sample);
      }
      sex = gr.getSex(sample);
      Diagnostic.userLog("Identified sex of sample " + sample + " from pedigree as " + sex);
    } else {
      sex = flags.isSet(MapFlags.SEX_FLAG) ? (Sex) flags.getValue(MapFlags.SEX_FLAG) : null;
    }
    return sex;
  }

  /**
   * Loads template and reads (left and right if paired end) in parallel.
   * @throws IOException if an IO problem occurs.
   */
  private static void makeSequenceParamsMulti(final NgsParamsBuilder ngsParamsBuilder, CFlags flags, final File build, final File buildSecond, NameParams nameParams, SequenceMode templateMode, SequenceMode readsMode) throws IOException {
    final LongRange buildRegion = CommonFlags.getReaderRestriction(flags);
    final ExecutorService executor = Executors.newFixedThreadPool(3);
    try {
      final FutureTask<SequenceParams> leftTask = new FutureTask<>(new SequenceParamsCallableSdf(build, buildRegion, nameParams, readsMode));
      executor.execute(leftTask);

      try {
        FutureTask<SequenceParams> rightTask = null;
        if (buildSecond != null) {
          rightTask = new FutureTask<>(new SequenceParamsCallableSdf(buildSecond, buildRegion, nameParams.includeFullNames() ? nameParams : new NameParams(false, false), readsMode));
          executor.execute(rightTask);
        }

        final FutureTask<SequenceParams> templateTask = getTemplateFutureTask(ngsParamsBuilder, flags, nameParams.includeFullNames(), templateMode);
        executor.execute(templateTask);

        ngsParamsBuilder.buildFirstParams(leftTask.get());
        if (rightTask != null) {
          ngsParamsBuilder.buildSecondParams(rightTask.get());
        }
        ngsParamsBuilder.searchParams(templateTask.get());
      } catch (final ExecutionException ie) {
        IOUtils.rethrow(ie.getCause());
      } catch (final InterruptedException ie) {
        throw new InvalidParamsException(ErrorType.INFO_ERROR, "Interrupted while loading datasets.");
      }
    } finally {
      executor.shutdownNow();
      while (!executor.isTerminated()) {
        try {
          executor.awaitTermination(1L, TimeUnit.SECONDS);
        } catch (InterruptedException ie) {
          break;
        }
      }
    }
  }

  /**
   * Loads template and reads (left and right if paired end) in parallel.
   * @throws IOException if an IO problem occurs.
   */
  private static void makeOnTheFlySequenceParamsMulti(NgsParamsBuilder ngsParamsBuilder, CFlags flags, DataSourceDescription desc, File build, File buildSecond, NameParams nameParams, boolean useQuality, SequenceMode templateMode, SequenceMode readsMode, SamSequenceReaderParams samParams) throws IOException {
    final FutureTask<SequenceParams> templateTask = getTemplateFutureTask(ngsParamsBuilder, flags, nameParams.includeFullNames(), templateMode);
    final LongRange buildRegion = CommonFlags.getReaderRestriction(flags);
    final ExecutorService executor = Executors.newFixedThreadPool(3);
    try {

      executor.execute(templateTask);

      try {
        final SimpleNames names = nameParams.includeNames() ? new SimpleNames() : null;
        final SimpleNames suffixes = nameParams.includeFullNames() ? new SimpleNames() : null;

        if (desc.getSourceFormat() == SourceFormat.SAM) {
          final FutureTask<SequenceParams[]> samTask = new FutureTask<>(new SequenceParamsCallableSam(build, desc, buildRegion, names, suffixes, useQuality, readsMode, samParams));
          executor.execute(samTask);
          final SequenceParams[] sp = samTask.get();
          assert sp.length == 2;
          ngsParamsBuilder.buildFirstParams(sp[0]);
          if (sp[1] != null) {
            ngsParamsBuilder.buildSecondParams(sp[1]);
          }
        } else {
          final FutureTask<SequenceParams> leftTask = new FutureTask<>(new SequenceParamsCallableFasta(build, desc, buildRegion, buildSecond != null ? PrereadArm.LEFT : PrereadArm.UNKNOWN, names, suffixes, useQuality, readsMode));
          executor.execute(leftTask);
          FutureTask<SequenceParams> rightTask = null;
          if (buildSecond != null) {
            final RightSimpleNames rNames = names == null ? null : new RightSimpleNames(names);
            final RightSimpleNames rSuffixes = suffixes == null ? null : new RightSimpleNames(suffixes);
            rightTask = new FutureTask<>(new SequenceParamsCallableFasta(buildSecond, desc, buildRegion, PrereadArm.RIGHT, rNames, rSuffixes, useQuality, readsMode));
            executor.execute(rightTask);
          }
          ngsParamsBuilder.buildFirstParams(leftTask.get());
          if (rightTask != null) {
            ngsParamsBuilder.buildSecondParams(rightTask.get());
          }
        }
        ngsParamsBuilder.searchParams(templateTask.get());
      } catch (final ExecutionException ie) {
        IOUtils.rethrow(ie.getCause());
      } catch (final InterruptedException ie) {
        throw new InvalidParamsException(ErrorType.INFO_ERROR, "Interrupted while loading datasets.");
      }
    } finally {
      executor.shutdownNow();
      while (!executor.isTerminated()) {
        try {
          executor.awaitTermination(1L, TimeUnit.SECONDS);
        } catch (InterruptedException ie) {
          break;
        }
      }
    }
  }

  /**
   * Retrieve the sam read group from somewhere within the flags object.
   * @param flags the flags object
   * @return the best sam read group record available from the flags
   * @throws IOException when reading one of the flag files fall
   */
  public static SAMReadGroupRecord getSAMReadGroupRecord(CFlags flags) throws IOException {
    final SAMReadGroupRecord rg;
    if (flags.isSet(SamCommandHelper.SAM_RG)) {
      final String readGroupString = (String) flags.getValue(SamCommandHelper.SAM_RG);
      rg = SamCommandHelper.validateAndCreateSamRG(readGroupString, SamCommandHelper.ReadGroupStrictness.REQUIRED);

    } else {
      final File reads = (File) flags.getValue(CommonFlags.READS_FLAG);
      final String formatName = (String) flags.getValue(FormatCli.FORMAT_FLAG);
      if (FormatCli.SDF_FORMAT.equals(formatName)) {
        final File sdf;
        if (ReaderUtils.isPairedEndDirectory(reads)) {
          sdf = ReaderUtils.getLeftEnd(reads);
        } else {
          sdf = reads;
        }
        final IndexFile index = new IndexFile(sdf);
        final String readGroupString = index.getSamReadGroup();
        if (readGroupString != null) {
          rg = SamCommandHelper.validateAndCreateSamRG(readGroupString.replaceAll("\t", "\\\\t"), SamCommandHelper.ReadGroupStrictness.REQUIRED);
        } else {
          rg = null;
        }
      } else if (SamCommandHelper.isSamInput(formatName) && reads.isFile()) {
        // Only try if it is a regular file (not a pipe which can't be opened twice)
        rg = SamCommandHelper.validateAndCreateSamRG(reads.getPath(), SamCommandHelper.ReadGroupStrictness.OPTIONAL);
      } else {
        rg = null;
      }
    }
    return rg;
  }

  static final class SequenceParamsCallableSdf implements Callable<SequenceParams> {
    private final File mBuild;
    private final boolean mUseMemReader;
    private final boolean mReads;
    private final LongRange mReaderRestriction;
    private final Sex mSex;
    private final boolean mIncludeReadNames;
    private final boolean mIncludeFullNames;
    private final SequenceMode mMode;

    SequenceParamsCallableSdf(final File build, final LongRange readerRestriction, NameParams nameParams, SequenceMode mode) { // C'tor for reads
      this(build, true, true, readerRestriction, null, nameParams.includeNames(), nameParams.includeFullNames(), mode);
    }

    SequenceParamsCallableSdf(File build, LongRange readerRestriction, boolean useMemReader, Sex sex, boolean includeFullNames, SequenceMode mode) { // C'tor for template
      this(build, useMemReader, false, readerRestriction, sex, false, includeFullNames, mode);
    }

    private SequenceParamsCallableSdf(File build, boolean useMemReader, boolean reads, LongRange readerRestriction, Sex sex, boolean includeReadNames, boolean includeFullNames, SequenceMode mode) {
      mBuild = build;
      mUseMemReader = useMemReader;
      mReads = reads;
      mReaderRestriction = readerRestriction;
      mSex = sex;
      mIncludeReadNames = includeReadNames;
      mIncludeFullNames = includeFullNames;
      mMode = mode;
    }

    @Override
    public SequenceParams call() {
      if (mReads) {
        return SequenceParams.builder().directory(mBuild).useMemReader(mUseMemReader).loadNames(mIncludeReadNames).loadFullNames(mIncludeFullNames).mode(mMode).readerRestriction(mReaderRestriction).create(); // Reads
      } else {
        SdfUtils.validateHasNames(mBuild);
        final SequenceParams params = SequenceParams.builder().sex(mSex).directory(mBuild).useMemReader(mUseMemReader).loadNames(true).loadFullNames(mIncludeFullNames).mode(mMode).readerRestriction(mReaderRestriction).create();
        try {
          SdfUtils.validateNoDuplicates(params.reader(), !mReads);
          return params; // Template
        } catch (final RuntimeException e) {
          try {
            params.close();
          } catch (final IOException e1) { }
          throw e;
        }
      }
    }
  }

  static final class SequenceParamsCallableFasta implements Callable<SequenceParams> {
    private final File mBuild;
    private final boolean mUseMemReader;
    private final DataSourceDescription mInputDescription;
    private final LongRange mReaderRestriction;
    private final SimpleNames mNames;
    private final SimpleNames mSuffixes;
    private final PrereadArm mArm;
    private final SequenceMode mMode;
    private final boolean mUseQuality;

    SequenceParamsCallableFasta(File build, DataSourceDescription desc, LongRange readerRestriction, PrereadArm arm, SimpleNames names, SimpleNames suffixes, boolean useQuality, SequenceMode mode) {
      // C'tor for reads
      mBuild = build;
      mUseMemReader = true;
      mReaderRestriction = readerRestriction;
      mInputDescription = desc;
      mNames = names;
      mSuffixes = suffixes;
      mArm = arm;
      mMode = mode;
      mUseQuality = useQuality;
    }

    @Override
    public SequenceParams call() throws IOException {
      final SequenceDataSource ds = FormatCli.getDnaDataSource(Collections.singletonList(mBuild), mInputDescription, mArm, false, false, null, false);
      final SequencesWriter sw = new SequencesWriter(ds, null, PrereadType.UNKNOWN, true);
      sw.setSdfId(new SdfId(0));
      final SequencesReader reader = sw.processSequencesInMemory(mBuild, mUseQuality, mNames, mSuffixes,  mReaderRestriction);
      return SequenceParams.builder().readerParam(new DefaultReaderParams(reader, mReaderRestriction, mMode)).useMemReader(mUseMemReader).mode(mMode).readerRestriction(mReaderRestriction).create(); // Reads
    }
  }

  static final class SequenceParamsCallableSam implements Callable<SequenceParams[]> {
    private final File mBuild;
    private final boolean mUseMemReader;
    private final DataSourceDescription mInputFormat;
    private final LongRange mReaderRestriction;
    private final SimpleNames mNames;
    private final SimpleNames mSuffixes;
    private final SequenceMode mMode;
    private final boolean mUseQuality;
    private final SamSequenceReaderParams mSamParams;

    SequenceParamsCallableSam(File build, DataSourceDescription desc, LongRange readerRestriction, SimpleNames names, SimpleNames suffixes, boolean useQuality, SequenceMode mode, SamSequenceReaderParams samParams) { // C'tor for reads
      mBuild = build;
      mUseMemReader = true;
      mReaderRestriction = readerRestriction;
      mInputFormat = desc;
      mNames = names;
      mSuffixes = suffixes;
      mMode = mode;
      mUseQuality = useQuality;
      mSamParams = samParams;
    }

    @Override
    public SequenceParams[] call() throws Exception {
      final SequenceDataSource ds = FormatCli.getDnaDataSource(Collections.singletonList(mBuild), mInputFormat, null, mSamParams.unorderedLoad(), mSamParams.flattenPairs(), null, false);
      final SequencesReader[] readers;
      if (mInputFormat.isInterleaved()) {
        final AlternatingSequencesWriter asw = new AlternatingSequencesWriter(ds, null, PrereadType.UNKNOWN, true);
        asw.setSdfId(new SdfId(0));
        asw.setCheckDuplicateNames(true);
        readers = asw.processSequencesInMemoryPaired(mBuild, mUseQuality, mNames, mSuffixes, mReaderRestriction);
      } else {
        final SequencesWriter sw = new SequencesWriter(ds, null, PrereadType.UNKNOWN, true);
        sw.setSdfId(new SdfId(0));
        sw.setCheckDuplicateNames(true);
        readers = new SequencesReader[] {sw.processSequencesInMemory(mBuild, mUseQuality, mNames, mSuffixes, mReaderRestriction), null};
      }
      final SequenceParams[] sp = new SequenceParams[readers.length];
      for (int i = 0; i < readers.length; ++i) {
        if (readers[i] != null) {
          sp[i] = SequenceParams.builder().readerParam(new DefaultReaderParams(readers[i], mReaderRestriction, mMode)).useMemReader(mUseMemReader).mode(mMode).readerRestriction(mReaderRestriction).create(); // Reads
        }
      }
      return sp;
    }
  }

  /**
   * Read aligner penalty flags into the {@link NgsParamsBuilder}
   * @param ngsParamsBuilder the builder
   * @param flags the flags containing settings
   * @return the builder for call chaining
   */
  public static NgsParamsBuilder populateAlignerPenaltiesParams(NgsParamsBuilder ngsParamsBuilder, CFlags flags) {
    ngsParamsBuilder.gapOpenPenalty((Integer) flags.getValue(MapFlags.GAP_OPEN_PENALTY_FLAG));
    ngsParamsBuilder.gapExtendPenalty((Integer) flags.getValue(MapFlags.GAP_EXTEND_PENALTY_FLAG));
    ngsParamsBuilder.substitutionPenalty((Integer) flags.getValue(MapFlags.MISMATCH_PENALTY_FLAG));
    ngsParamsBuilder.unknownsPenalty((Integer) flags.getValue(MapFlags.UNKNOWNS_PENALTY_FLAG));
    ngsParamsBuilder.indelSoftClipDistance((Integer) flags.getValue(MapFlags.SOFT_CLIP_DISTANCE_FLAG));
    ngsParamsBuilder.alignerBandWidthFactor(new MaxShiftFactor((Double) flags.getValue(MapFlags.ALIGNER_BAND_WIDTH_FACTOR_FLAG)));

    ngsParamsBuilder.mismatchSoftClipDistance(GlobalFlags.getIntegerValue(CoreGlobalFlags.EDIT_DIST_MISMATCH_SOFT_CLIP));
    ngsParamsBuilder.minMatches(GlobalFlags.getIntegerValue(CoreGlobalFlags.EDIT_DIST_MIN_MATCHES));
    return ngsParamsBuilder;
  }
}
