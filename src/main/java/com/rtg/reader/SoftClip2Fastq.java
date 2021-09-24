/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.PAIR_ORIENTATION_FLAG;
import static com.rtg.reader.Sdf2Fasta.INTERLEAVE;
import static com.rtg.reader.Sdf2Fasta.LINE_LENGTH;
import static com.rtg.reader.Sdf2Fasta.OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;
import static com.rtg.util.machine.MachineOrientation.FR;
import static com.rtg.util.machine.MachineOrientation.RF;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNA;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamFilterOptions;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.BaseFile;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineOrientation;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Finds reads with long soft clips, and splits the read at that point into two parts for re-mapping.
 */
public class SoftClip2Fastq extends AbstractCli {

  private static final String MIN_SOFT_CLIP_LENGTH = "min-soft-clip-length";
  private static final String SOFT_CLIPS_ONLY = "soft-clipped-only";
  private static final String EXCLUDE_SECONDARY = "exclude-secondary";
  private static final String REVCOMP = "reverse-complement-r2";

  private int mSplitRecords;
  private MachineOrientation mOrientation;

  @Override
  public String moduleName() {
    return "softclip2fastq";
  }

  @Override
  public String description() {
    return "extract long soft clips from read alignments as FASTQ";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Select alignments containing long soft clips and output as FASTQ. By default the output is to create a pseudo read pair, the first read corresponds to the soft clipped bases and the second to the aligned bases.");
    CommonFlagCategories.setCategories(mFlags);

    final Flag<File> inFlag = mFlags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing coordinate-sorted alignments")
      .setCategory(INPUT_OUTPUT).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing coordinate-sorted alignments").setCategory(INPUT_OUTPUT);
    SamFilterOptions.registerRestrictionFlag(mFlags);
    SamFilterOptions.registerBedRestrictionFlag(mFlags);
    Sdf2Fasta.registerTextOutputFlags(mFlags);
    mFlags.registerOptional(EXCLUDE_SECONDARY, "if set, exclude secondary alignments").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_SOFT_CLIP_LENGTH, Integer.class, CommonFlags.INT, "minimum length of soft clip to trigger read splitting", 35).setCategory(SENSITIVITY_TUNING);

    mFlags.registerOptional(SOFT_CLIPS_ONLY, "output only the portion of the read corresponding to the soft-clipped bases. Default is to output as a pseudo-read-pair").setCategory(UTILITY);
    mFlags.registerOptional(INTERLEAVE, "interleave paired data into a single output file. Default is to split to separate output files").setCategory(UTILITY);
    mFlags.registerOptional(PAIR_ORIENTATION_FLAG, MachineOrientation.class, CommonFlags.STRING, "orientation for pseudo pairs", FR).setCategory(UTILITY);

    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);

    mFlags.setValidator(flags -> {
        final File baseOutput = (File) flags.getValue(OUTPUT);
        final boolean gzip = !mFlags.isSet(NO_GZIP);
        final BaseFile baseFile = FastqUtils.baseFile(baseOutput, gzip);
        return CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
          && (flags.isSet(SOFT_CLIPS_ONLY) || PairedEndTrimCli.validatedPairedOutputFiles(flags, baseOutput, baseFile))
          && SamFilterOptions.validateFilterFlags(flags, true)
          && Sdf2Fasta.validateTextOutputFlags(flags)
          && flags.checkInRange(MIN_SOFT_CLIP_LENGTH, 1, Integer.MAX_VALUE)
          ;
      }
    );
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream log) throws IOException {
    final File baseOutput = (File) mFlags.getValue(OUTPUT);
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final int lineLength = (Integer) mFlags.getValue(LINE_LENGTH);
    final int minSoftClipLength = (Integer) mFlags.getValue(MIN_SOFT_CLIP_LENGTH);
    final boolean softClipSideOnly = mFlags.isSet(SOFT_CLIPS_ONLY);
    final boolean interleavePaired = mFlags.isSet(INTERLEAVE);
    if (!interleavePaired && !softClipSideOnly && FileUtils.isStdio(baseOutput)) {
      throw new NoTalkbackSlimException("Sending non-interleaved paired-end data to stdout is not supported.");
    }
    mOrientation = (MachineOrientation) mFlags.getValue(PAIR_ORIENTATION_FLAG);

    final SamFilterParams.SamFilterParamsBuilder fpb = SamFilterOptions.makeFilterParamsBuilder(mFlags);
    if (mFlags.isSet(EXCLUDE_SECONDARY)) {
      fpb.requireUnsetFlags(SamBamConstants.SAM_SECONDARY_ALIGNMENT);
    }
    final Collection<File> inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS, CommandLineFiles.NOT_DIRECTORY).getFileList(mFlags);
    final int threads = CommonFlags.parseIOThreads(null); // Just use the default
    final SamReadingContext c = new SamReadingContext(inputFiles, threads, fpb.create(), SamUtils.getUberHeader(inputFiles), null);
    try (RecordIterator<SAMRecord> reader = new ThreadedMultifileIterator<>(c, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      final BaseFile baseFile = FastqUtils.baseFile(baseOutput, gzip);
      try (final FastqWriter left = interleavePaired || softClipSideOnly
        ? new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "")), lineLength, (byte) 0, false)
        : new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "_1")), lineLength, (byte) 0, false);
           final FastqWriter right = softClipSideOnly ? null : interleavePaired
             ? left
             : new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "_2")), lineLength, (byte) 0, false)) {
        while (reader.hasNext()) {
          if (checkRecord(reader.next(), minSoftClipLength, left, right)) {
            ++mSplitRecords;
          }
        }

      }
    }
    if (!FileUtils.isStdio(baseOutput)) {
      Diagnostic.info("Split records: " + mSplitRecords);
    }
    return 0;
  }

  private boolean checkRecord(SAMRecord record, int minSoftClipLength, FastqWriter left, FastqWriter right) throws IOException {
    final List<CigarElement> cigar = record.getCigar().getCigarElements();
    if (cigar.size() > 1) {
      final CigarElement first = cigar.get(0);
      if ((first.getOperator() == CigarOperator.SOFT_CLIP) && (first.getLength() >= minSoftClipLength)) {
        splitAndWrite(left, right, record, first.getLength(), true);
        return true;
      } else {
        final CigarElement last = cigar.get(cigar.size() - 1);
        if ((last.getOperator() == CigarOperator.SOFT_CLIP) && (last.getLength() >= minSoftClipLength)) {
          splitAndWrite(left, right, record, record.getReadLength() - last.getLength(), false);
          return true;
        }
      }
    }
    return false;
  }

  private void splitAndWrite(FastqWriter left, FastqWriter right, SAMRecord record, int splitPos, boolean clippedOnStart) throws IOException {
    final byte[] seq = DNA.byteDNAtoByte(record.getReadBases());
    final byte[] qual = record.getBaseQualities();
    final String info = (record.getReadPairedFlag() ? record.getFirstOfPairFlag() ? "_1 " : "_2 " : " ")
      + record.getReferenceName() + ":" + record.getAlignmentStart() + "/" + record.getFlags() + "/" + record.getCigarString();
    FastqSequence l = new FastqSequence(record.getReadName() + info + "/start", seq, qual, splitPos);
    FastqSequence r = new FastqSequence(record.getReadName() + info + "/end", Arrays.copyOfRange(seq, splitPos, seq.length), Arrays.copyOfRange(qual, splitPos, seq.length), seq.length - splitPos);
    if (!clippedOnStart) { // Swap to ensure that l corresponds to the soft clipped sequence
      final FastqSequence t = l; l = r; r = t;
      // We should rc both arms here, but we'll do it below to avoid the work of doing a double-rc (depending on orientation)
    }

    if (!clippedOnStart ^ mOrientation == RF) {  // Flip first side if we swapped or if we want the resulting pseudo-pair should align RF style (but not both)
      l.rc();
    }
    left.write(l);
    if (right != null) {
      if (!clippedOnStart ^ mOrientation == FR) {  // Flip second side if we swapped or if we want the resulting pseudo-pair should align FR style (but not both)
        r.rc();
      }
      right.write(r);
    }
  }
}
