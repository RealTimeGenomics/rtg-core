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

package com.rtg.variant.sv;

import static com.rtg.launcher.CommonFlags.FORCE;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.BamIndexer;
import com.rtg.sam.SamUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

/**
 * Augment unmated SAM files and output read group statistics
 */
public class UnmatedAugmenterCli extends LoggedCli {

  static final String KEEP_ORIG_FLAG = "Xkeep-original";
  static final String OUTPUT_SUFFIX_FLAG = "Xoutput-suffix";
  static final String NO_AUGMENT = "no-augment";

  @Override
  public String moduleName() {
    return "svprep";
  }

  @Override
  public String description() {
    return "prepare SAM/BAM files for structural variant analysis";
  }

  @Override
  protected void mainExit(String[] args) {
    super.mainExit(args);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getAnonymousValue(0);
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Prepares SAM files for use with sv module, by finding mate position/frame information for discordant matings and generates the read group statistics. This command is not needed with RTG mappings unless automatic svprep was disabled.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerExtendedHelp();
    CommonFlags.initForce(mFlags);
    mFlags.registerOptional('s', OUTPUT_SUFFIX_FLAG, String.class, CommonFlags.STRING, "suffix for output file of each input file", ".augmented").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('k', KEEP_ORIG_FLAG, "keep original file and create augmented file using suffix").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(NO_AUGMENT, "if set, only compute read group statistics").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerRequired(File.class, CommonFlags.DIR, "directory containing SAM/BAM format files").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.setValidator(flags -> checkInputDirectory(flags) && flags.checkIf(OUTPUT_SUFFIX_FLAG, KEEP_ORIG_FLAG));
  }

  private static boolean checkInputDirectory(CFlags flags) {
    final File dir = (File) flags.getAnonymousValue(0);
    if (!dir.isDirectory()) {
      flags.setParseMessage("A directory name was expected, but \"" + dir.getPath() + "\" is not a directory.");
      return false;
    }
    return true;
  }

  static class AlignmentsFileFilter implements FileFilter {
    @Override
    public boolean accept(File pathname) {
      final String fname = pathname.getName().toLowerCase(Locale.ROOT);
      return fname.matches(".*\\.sam(.gz)?") || fname.matches(".*\\.bam");
    }
  }
  static class MatedFileFilter implements FileFilter {
    @Override
    public boolean accept(File pathname) {
      final String fname = pathname.getName().toLowerCase(Locale.ROOT);
      return fname.matches("mated\\.sam(.gz)?") || fname.matches("mated\\.bam");
    }
  }
  static class UnmatedFileFilter implements FileFilter {
    @Override
    public boolean accept(File pathname) {
      final String fname = pathname.getName().toLowerCase(Locale.ROOT);
      return fname.matches("unmated\\.sam(.gz)?") || fname.matches("unmated\\.bam");
    }
  }
  static class UnmappedFileFilter implements FileFilter {
    @Override
    public boolean accept(File pathname) {
      final String fname = pathname.getName().toLowerCase(Locale.ROOT);
      return fname.matches("unmapped\\.sam(.gz)?") || fname.matches("unmapped\\.bam");
    }
  }

  private File getFile(File dir, FileFilter filter, String label) throws IOException {
    final File[] files = FileUtils.listFiles(dir, filter);
    if (files.length > 1) {
      throw new NoTalkbackSlimException("Multiple " + label + " files found!");
    }
    return files.length == 0 ? null : files[0];
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File dir = outputDirectory();
    final File matedFile = getFile(dir, new MatedFileFilter(), "mated");
    final File unmatedFile = getFile(dir, new UnmatedFileFilter(), "unmated");
    final File unmappedFile = getFile(dir, new UnmappedFileFilter(), "unmapped");
    final File[] alignmentsFiles = FileUtils.listFiles(dir, new AlignmentsFileFilter());
    if (matedFile == null && unmatedFile == null && alignmentsFiles.length == 0) {
      throw new NoTalkbackSlimException("Could not find any alignments files to process");
    }
    if (matedFile == null ^ unmatedFile == null) {
      throw new NoTalkbackSlimException("Could not find both mated and unmated alignments files.");
    }

    try (OutputStream rgOut = FileUtils.createOutputStream(new File(dir, UnmatedAugmenter.DEFAULT_RGSTATS_FILENAME))) {
      if (mFlags.isSet(NO_AUGMENT)) {
        new ReadGroupStatsCalculator().calculate(Arrays.asList(alignmentsFiles), rgOut);
      } else if (matedFile != null && unmatedFile != null) {
        augmentSplitFiles(matedFile, unmatedFile, unmappedFile, !mFlags.isSet(KEEP_ORIG_FLAG), mFlags.isSet(FORCE), (String) mFlags.getValue(OUTPUT_SUFFIX_FLAG), rgOut);
      } else {
        augmentMergedFiles(alignmentsFiles, !mFlags.isSet(KEEP_ORIG_FLAG), mFlags.isSet(FORCE), (String) mFlags.getValue(OUTPUT_SUFFIX_FLAG), rgOut);
      }
    }
    return 0;
  }

  /**
   * Augments self-contained alignment files with potential mate information
   * @param files files containing mated/unmated/unmapped alignments
   * @param rename rename result file over top of original
   * @param force overwrite any files that match output file name derived from input and suffix. Ignored when rename set to true.
   * @param suffix suffix to add to filename appears before optional <code>.sam[.gz]</code>.
   * @param readGroupStatsOut stream to write out read group stats
   * @throws IOException they happen maybe
   */
  static void augmentMergedFiles(File[] files, boolean rename, boolean force, String suffix, OutputStream readGroupStatsOut) throws IOException {
    final ReadGroupStatsCalculator outer = new ReadGroupStatsCalculator();

    for (File file : files) {
      final ReadGroupStatsCalculator calc = new ReadGroupStatsCalculator();
      final UnmatedAugmenter augmenter = new UnmatedAugmenter();
      final File output = getOutputFile(file, rename, force, suffix);
      augmenter.augmentMixed(file, output, calc);
      renameAndReindex(file, output, rename);
      Diagnostic.userLog("Agumented " + augmenter.mAugmentedUnmated + " unmated and " + augmenter.mAugmentedUnmapped + " unmapped records");
      outer.merge(calc);
    }

    outer.writeReadGroupStats(readGroupStatsOut);
  }


  /**
   * Augments the files with potential mate information
   * @param matedFile mated file
   * @param unmatedFile unmated file
   * @param unmappedFile unmapped file
   * @param rename rename result file over top of original
   * @param force overwrite any files that match output file name derived from input and suffix. Ignored when rename set to true.
   * @param suffix suffix to add to filename appears before optional <code>.sam[.gz]</code>.
   * @param readGroupStatsOut stream to write out read group stats
   * @throws IOException they happen maybe
   */
  static void augmentSplitFiles(File matedFile, File unmatedFile, File unmappedFile, boolean rename, boolean force, String suffix, OutputStream readGroupStatsOut) throws IOException {
    final ReadGroupStatsCalculator calc = new ReadGroupStatsCalculator();
    calc.addFile(matedFile);

    final File output = getOutputFile(unmatedFile, rename, force, suffix);
    final UnmatedAugmenter augmenter = new UnmatedAugmenter();
    augmenter.augmentUnmated(unmatedFile, output, calc);
    renameAndReindex(unmatedFile, output, rename);

    if (unmappedFile != null) {
      final File outputUnmapped = getOutputFile(unmappedFile, rename, force, suffix);
      augmenter.augmentUnmapped(unmappedFile, outputUnmapped, calc);
      renameAndReindex(unmappedFile, outputUnmapped, rename);
    }

    Diagnostic.userLog("Agumented " + augmenter.mAugmentedUnmated + " unmated and " + augmenter.mAugmentedUnmapped + " unmapped records");
    calc.writeReadGroupStats(readGroupStatsOut);
  }

  private static void renameAndReindex(File dest, File src, boolean rename) throws IOException {
    final File finalFile;
    if (rename) {
      if (!(dest.canWrite() && dest.delete() && src.renameTo(dest))) {
        throw new NoTalkbackSlimException("Unable to rename file: \"" + src.getPath() + "\" to: \"" + dest.getPath() + "\"");
      }
      finalFile = dest;
    } else {
      finalFile = src;
    }

    //recreate index if it exists
    File indexFile;
    if (SamUtils.isBAMFile(finalFile)) {
      indexFile = BamIndexer.indexFileName(finalFile);
      if (!indexFile.exists()) {
        indexFile = BamIndexer.secondaryIndexFileName(finalFile);
      }
      if (indexFile.exists()) {
        try {
          BamIndexer.saveBamIndex(finalFile, indexFile);
        } catch (final UnindexableDataException e) {
          Diagnostic.warning("Could not create BAM index: " + indexFile.getPath());
          if (indexFile.exists() && !indexFile.delete()) {
            throw new IOException("Could not delete old index file: " + indexFile.getPath());
          }
        }
      }
    } else {
      indexFile = TabixIndexer.indexFileName(finalFile);
      if (indexFile.exists()) {
        try {
          new TabixIndexer(finalFile, indexFile).saveSamIndex();
        } catch (final UnindexableDataException e) {
          Diagnostic.warning("Could not create SAM index: " + indexFile.getPath());
          if (indexFile.exists() && !indexFile.delete()) {
            throw new IOException("Could not delete old index file: " + indexFile.getPath());
          }
        }
      }
    }
  }

  static File getOutputFile(File f, boolean rename, boolean force, String suffix) throws IOException {
    final String name = f.getName().replaceAll("^(.*?)(\\.[bs]am(\\.gz)?)?$", "$1" + suffix + "$2");
    File ret = new File(f.getParentFile(), name);
    if (rename && ret.exists()) {
      ret = File.createTempFile(f.getName(), ".temp", f.getParentFile());
    } else if (ret.exists() && !force) {
      throw new NoTalkbackSlimException("File: \"" + ret.getPath() + "\" already exists. Please remove or choose a different suffix.");
    }
    return ret;
  }

}
