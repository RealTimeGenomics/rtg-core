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

package com.rtg.visualization;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.FlagValue;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.intervals.RegionRestriction;


/**
 * Hold parameters for aview
 */
@TestClass(value = {"com.rtg.visualization.AviewParamsTest", "com.rtg.visualization.AviewTest"})
class AviewParams {

  static final String BASELINE_VCF = "baseline";
  static final String REFERENCE = "template";
  static final String CALLS_VCF = "calls";
  static final String BEDS = "bed";
  static final String ALIGNMENTS = "alignments";

  static final String SAMPLE = "sample";
  static final String REGION = "region";
  static final String MAX_MATED_ALIGNMENT_SCORE = "Xmax-as-mated";
  static final String MAX_UNMATED_ALIGNMENT_SCORE = "Xmax-as-unmated";
  static final String MAX_IH_SCORE = "Xmax-hits";
  static final String MIN_MAPQ_SCORE = "Xmin-mapq";

  static final String PRINT_CIGARS = "print-cigars";
  static final String PRINT_READGROUP = "print-readgroup";
  static final String PRINT_MAPQ = "print-mapq";
  static final String PRINT_NAMES = "print-names";
  static final String PRINT_MATE_POSITION = "print-mate-position";
  static final String PRINT_REFERENCE_LINE = "print-reference-line";

  static final String NO_DOTS = "no-dots";
  static final String NO_COLOR = "no-color";
  static final String HTML = "html";
  static final String PADDING = "padding";
  static final String SORT_READGROUP = "sort-readgroup";
  static final String SORT_READS = "sort-reads";
  static final String PROJECT_TRACK = "project-track";
  static final String UNFLATTEN = "unflatten";

  static final String READS_SDF = "reads";
  static final String UNMAPPED_SAM = "Xunmapped";
  static final String XMAPPING_TOLERANCE = "Xmapping-tolerance";


  static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);

    // Input sources
    flags.registerRequired('t', REFERENCE, File.class, "SDF", "reference SDF").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('r', READS_SDF, File.class, "SDF", "read SDF (only needed to indicate correctness of simulated read mappings)")
        .setMaxCount(Integer.MAX_VALUE).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag in = flags.registerRequired(File.class, "FILE", "alignment SAM/BAM files").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    in.setMinCount(0);
    in.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line)").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('b', BASELINE_VCF, File.class, "FILE", "VCF file containing baseline variants").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('c', CALLS_VCF, File.class, "FILE", "VCF file containing called variants").setMaxCount(Integer.MAX_VALUE).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('B', BEDS, File.class, "FILE", "BED file containing regions to overlay").setMaxCount(Integer.MAX_VALUE).setCategory(CommonFlagCategories.INPUT_OUTPUT);

    // Filtering options
    flags.registerOptional(SAMPLE, String.class, "STRING", "specify name of sample to select").setMaxCount(Integer.MAX_VALUE).setCategory(CommonFlagCategories.FILTERING);
    flags.registerRequired(REGION, String.class, "STRING", "the region of interest to display. The format is one of <sequence_name>, <sequence_name>:start-end or <sequence_name>:start+length").setCategory(CommonFlagCategories.FILTERING);
    flags.registerOptional('p', PADDING, Integer.class, "INT", "padding around region of interest. Default is to automatically determine padding to avoid read truncation").setCategory(CommonFlagCategories.FILTERING);

    flags.registerOptional(MAX_MATED_ALIGNMENT_SCORE, Integer.class, "INT", "if set, ignore mated SAM records with an alignment score (AS attribute) that exceeds this value").setCategory(CommonFlagCategories.FILTERING);
    flags.registerOptional(MAX_UNMATED_ALIGNMENT_SCORE, Integer.class, "INT", "if set, ignore unmated SAM records with an alignment score (AS attribute) that exceeds this value").setCategory(CommonFlagCategories.FILTERING);
    flags.registerOptional(MAX_IH_SCORE, Integer.class, "INT", "if set, ignore SAM records with a number of hits (NH attribute) that exceeds this value").setCategory(CommonFlagCategories.FILTERING);
    flags.registerOptional(MIN_MAPQ_SCORE, Integer.class, "INT", "if set, ignore SAM records with a MAPQ score less than this value").setCategory(CommonFlagCategories.FILTERING);

    // What to display and how
    flags.registerOptional(PRINT_CIGARS, "print alignment cigars").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(PRINT_READGROUP, "print read group id for each alignment").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(PRINT_MATE_POSITION, "print mate position").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(PRINT_MAPQ, "print alignment MAPQ values").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(PRINT_NAMES, "print read names").setCategory(CommonFlagCategories.REPORTING);

    // Other display options
    flags.registerOptional(NO_DOTS, "display nucleotide instead of dots").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(NO_COLOR, "do not use colors").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(HTML, "output as HTML").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(PRINT_REFERENCE_LINE, Integer.class, "INT", "print reference line every N lines", 0).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(PROJECT_TRACK, Integer.class, "INT", "if set, project highlighting for the specified track down through reads. Default projects the union of tracks").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(SORT_READS, "sort reads on start position").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(SORT_READGROUP, "sort reads first on read group and then on start position").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(UNFLATTEN, "display unflattened CGI reads when present").setCategory(CommonFlagCategories.REPORTING);

    // Random extras
    flags.registerOptional(XMAPPING_TOLERANCE, Integer.class, "INT", "variation allowed in start position when determining correctness of simulated read mapping", 0).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional('U', UNMAPPED_SAM, File.class, "FILE", "unmapped SAM file (only for simulated data)")
        .setMaxCount(Integer.MAX_VALUE).setCategory(CommonFlagCategories.REPORTING);

    flags.addRequiredSet(in);
    flags.addRequiredSet(listFlag);
    flags.setValidator(new Validator() {

      @Override
      public boolean isValid(CFlags flags) {
        // It's OK to not provide mappings, but we need at least one type of file to display
        if (!(CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE) || flags.isSet(BASELINE_VCF) || flags.isSet(CALLS_VCF) || flags.isSet(BEDS))) {
          return false;
        }
        if (flags.isSet(SORT_READS) && flags.isSet(SORT_READGROUP)) {
          flags.setParseMessage("only one of --" + SORT_READS + ", or --" + SORT_READGROUP + " must be set.");
          return false;
        }
        final String region = (String) flags.getValue(REGION);
        try {
          new RegionRestriction(region);
        } catch (IllegalArgumentException e) {
          flags.setParseMessage("invalid --" + REGION + " specification " + region);
          return false;
        }
        final int variance = (Integer) flags.getValue(XMAPPING_TOLERANCE);
        if (variance < 0) {
          Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + XMAPPING_TOLERANCE, variance + "", "0");
          return false;
        }
        return true;
      }
    });
  }

  private final File[] mAlignments;
  private final File[] mUnmapped;
  private final File[] mTrackFiles;
  private final File mBaselineFile;
  private final File mReference;
  private final File[] mReads;
  private final String[] mSamples;
  private final boolean mDisplayDots;
  private final boolean mUseTerminalColor;
  private final boolean mUseHtml;

  private final RegionRestriction mRegion;
  private final int mRegionPadding;
  private final int mHeaderLineRepeat;
  private final int mMaxMatedAlignmentScore;
  private final int mMaxUnmatedAlignmentScore;
  private final int mMaxIhScore;
  private final int mMinMapq;
  private final int mMappingTolerance;
  private final boolean mPrintCigars;
  private final boolean mPrintReadGroup;
  private final boolean mPrintMatePosition;
  private final boolean mPrintReadName;
  private final boolean mSortReads;
  private final boolean mSortReadGroup;
  private final boolean mPrintMapQ;
  private final boolean mUnflattenCgi;
  private final int mProjectTrackId;

  AviewParams(AviewParamsBuilder aviewParamsBuilder) {
    mAlignments = aviewParamsBuilder.mAlignments;
    mDisplayDots = aviewParamsBuilder.mDisplayDots;
    // We will treat a region of foobar:start the same as foobar:start+1
    mRegion = aviewParamsBuilder.mRegion.getEnd() == RegionRestriction.MISSING && aviewParamsBuilder.mRegion.getStart() != 0
        ? new RegionRestriction(aviewParamsBuilder.mRegion.getSequenceName(), aviewParamsBuilder.mRegion.getStart(), aviewParamsBuilder.mRegion.getStart() + 1)
        : aviewParamsBuilder.mRegion;
     mRegionPadding = aviewParamsBuilder.mRegionPadding;
    mBaselineFile = aviewParamsBuilder.mBaselineFile;
    mTrackFiles = aviewParamsBuilder.mTrackFiles;
    mHeaderLineRepeat = aviewParamsBuilder.mHeaderLineRepeat;
    mMaxMatedAlignmentScore = aviewParamsBuilder.mMaxMatedAlignmentScore;
    mMaxUnmatedAlignmentScore = aviewParamsBuilder.mMaxUnmatedAlignmentScore;
    mMaxIhScore = aviewParamsBuilder.mMaxIhScore;
    mMinMapq = aviewParamsBuilder.mMinMapq;
    mMappingTolerance = aviewParamsBuilder.mMappingTolerance;
    mUseTerminalColor = aviewParamsBuilder.mUseTerminalColor;
    mUseHtml = aviewParamsBuilder.mUseHtml;
    mPrintCigars = aviewParamsBuilder.mPrintCigars;
    mPrintReadName = aviewParamsBuilder.mPrintReadName;
    mReads = aviewParamsBuilder.mReads;
    mReference = aviewParamsBuilder.mReference;
    mUnmapped = aviewParamsBuilder.mUnmapped;
    mSortReads = aviewParamsBuilder.mSortReads;
    mSortReadGroup = aviewParamsBuilder.mSortReadGroup;
    mPrintReadGroup = aviewParamsBuilder.mPrintReadGroup;
    mPrintMatePosition = aviewParamsBuilder.mPrintMatePosition;
    mSamples = aviewParamsBuilder.mSamples;
    mPrintMapQ = aviewParamsBuilder.mPrintMapQ;
    mUnflattenCgi = aviewParamsBuilder.mUnflattenCgi;
    mProjectTrackId = aviewParamsBuilder.mProjectTrackId;
  }

 private static File[] getFiles(CFlags flags, String flagName) {
    final ArrayList<File> list = new ArrayList<>();
    final Collection<Object> files;
    if (ALIGNMENTS.equals(flagName)) {
      files = flags.getAnonymousValues(0);
    } else {
      files = flags.getValues(flagName);
    }
    for (final Object file : files) {
      list.add((File) file);
    }
    return list.toArray(new File[list.size()]);
  }

 static AviewParams makeParams(CFlags flags) throws IOException {
   final AviewParamsBuilder builder = new AviewParamsBuilder();
   final List<File> fileList = CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, false);
   File[] list = new File[fileList.size()];
   list = fileList.toArray(list);
   builder.alignments(list)
     .reference((File) flags.getValue(REFERENCE))
     .region(flags.getValue(REGION).toString())
     .headerLineRepeat((Integer) flags.getValue(PRINT_REFERENCE_LINE))
     .displayDots(!flags.isSet(NO_DOTS))
     .useTerminalColor(!flags.isSet(NO_COLOR))
     .printCigars(flags.isSet(PRINT_CIGARS))
     .printReadName(flags.isSet(PRINT_NAMES))
     .printMapQ(flags.isSet(PRINT_MAPQ))
     .useHtml(flags.isSet(HTML))
     .unflattenCgi(flags.isSet(UNFLATTEN))
     .sortReads(flags.isSet(SORT_READS))
     .sortReadGroup(flags.isSet(SORT_READGROUP))
     .printReadGroup(flags.isSet(PRINT_READGROUP))
     .printMatePosition(flags.isSet(PRINT_MATE_POSITION))
     ;
   if (flags.isSet(MAX_IH_SCORE)) {
     builder.maxIhScore((Integer) flags.getValue(MAX_IH_SCORE));
   }
   if (flags.isSet(MIN_MAPQ_SCORE)) {
     builder.minMapq((Integer) flags.getValue(MIN_MAPQ_SCORE));
   }
   if (flags.isSet(MAX_MATED_ALIGNMENT_SCORE)) {
     builder.maxMatedAlignmentScore((Integer) flags.getValue(MAX_MATED_ALIGNMENT_SCORE));
   }
   if (flags.isSet(MAX_UNMATED_ALIGNMENT_SCORE)) {
     builder.maxUnmatedAlignmentScore((Integer) flags.getValue(MAX_UNMATED_ALIGNMENT_SCORE));
   }
   if (flags.isSet(READS_SDF)) {
     builder.reads(getFiles(flags, READS_SDF));
   }
   if (flags.isSet(UNMAPPED_SAM)) {
     builder.unmapped(getFiles(flags, UNMAPPED_SAM));
   }
   if (flags.isSet(PADDING)) {
     builder.regionPadding((Integer) flags.getValue(PADDING));
   }
   if (flags.isSet(PROJECT_TRACK)) {
     builder.projectTrackId((Integer) flags.getValue(PROJECT_TRACK));
   }

   if (flags.isSet(BASELINE_VCF)) {
     builder.baselineFile((File) flags.getValue(BASELINE_VCF));
   }

   final ArrayList<File> tracks = new ArrayList<>();
   final Flag callFlag = flags.getFlag(CALLS_VCF);
   final Flag bedFlag = flags.getFlag(BEDS);
   for (final FlagValue flag : flags.getReceivedValues()) {
     if (flag.getFlag() == bedFlag || flag.getFlag() == callFlag) {
       tracks.add((File) flag.getValue());
     }
   }
   builder.trackFiles(tracks.toArray(new File[tracks.size()]));

   if (flags.isSet(SAMPLE)) {
     final ArrayList<String> s = new ArrayList<>();
     for (Object o : flags.getValues(SAMPLE)) {
       s.add((String) o);
     }
     builder.samples(s.toArray(new String[s.size()]));
   }
   builder.mappingTolerance((Integer) flags.getValue(XMAPPING_TOLERANCE));

   return builder.create();
 }

  final File[] readFiles() {
    return mReads;
  }

  final File[] alignmentsFiles() {
    return mAlignments;
  }

  final File[] unmappedFiles() {
    return mUnmapped;
  }

  final File[] trackFiles() {
    return mTrackFiles;
  }

  final File baselineFile() {
    return mBaselineFile;
  }

  final String[] wantedSamples() {
    return mSamples;
  }

  final File referenceFile() {
    return mReference;
  }


  final String sequenceName() {
    return mRegion.getSequenceName();
  }

  final RegionRestriction region() {
    return mRegion;
  }

  // one based inclusive
  final int start() {
    return mRegion.getStart() == RegionRestriction.MISSING ? -1 : mRegion.getStart() + 1;
  }

  // Now one based exclusive
  final int end() {
    return mRegion.getStart() == RegionRestriction.MISSING ? -1 : mRegion.getEnd() + 1;
  }

  final int regionPadding() {
    return mRegionPadding;
  }

  final int headerLineRepeat() {
    return mHeaderLineRepeat;
  }


  final int maxMatedAlignmentScore() {
    return mMaxMatedAlignmentScore;
  }

  final int maxUnmatedAlignmentScore() {
    return mMaxUnmatedAlignmentScore;
  }

  final int minMapq() {
    return mMinMapq;
  }

  final int maxIhScore() {
    return mMaxIhScore;
  }

  final boolean displayDots() {
    return mDisplayDots;
  }

  final boolean printCigars() {
    return mPrintCigars;
  }

  final boolean printReadName() {
    return mPrintReadName;
  }

  final boolean useTerminalColor() {
    return mUseTerminalColor;
  }

  final boolean unflattenCgi() {
    return mUnflattenCgi;
  }

  final int projectTrackId() {
    return mProjectTrackId;
  }

  final DisplayHelper displayHelper() {
    return mUseHtml
    ? new HtmlDisplayHelper()
    : mUseTerminalColor
    ? new AnsiDisplayHelper()
    : new DisplayHelper();
  }

  final boolean sortReads() {
    return mSortReads;
  }

  final boolean sortReadGroup() {
    return mSortReadGroup;
  }

  final boolean printReadGroup() {
    return mPrintReadGroup;
  }

  final boolean printMapQ() {
    return mPrintMapQ;
  }

  final boolean printMatePosition() {
    return mPrintMatePosition;
  }

  final int mappingTolerance() {
    return mMappingTolerance;
  }
}
