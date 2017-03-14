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

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.BamIndexer;
import com.rtg.sam.SamUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
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
    mFlags.registerOptional('s', OUTPUT_SUFFIX_FLAG, String.class, CommonFlags.STRING, "suffix for output file of each input file", ".augmented").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('k', KEEP_ORIG_FLAG, "keep original file and create augmented file using suffix").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired(File.class, CommonFlags.DIR, "directory containing SAM/BAM format files").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.setValidator(new UnmatedAugmenterValidator());
  }

  static class SamBamGzMatedUnmatedFileFilter implements FileFilter {
    @Override
    public boolean accept(File pathname) {
      final String fname = pathname.getName().toLowerCase(Locale.ROOT);
      return fname.matches(".*?mated\\.sam(.gz)?") || fname.matches(".*?mated\\.bam");
    }
  }

  static class SamBamGzUnmappedFileFilter implements FileFilter {
    @Override
    public boolean accept(File pathname) {
      final String fname = pathname.getName().toLowerCase(Locale.ROOT);
      return fname.matches(".*?unmapped\\.sam(.gz)?") || fname.matches(".*?unmapped\\.bam");
    }
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File dir = outputDirectory();
    final File[] files = dir.listFiles(new SamBamGzMatedUnmatedFileFilter());
    if (files == null) {
      throw new IOException("Invalid directory.");
    } else if (files.length == 0) {
      Diagnostic.warning("Read group statistics calculation requires properly mated data");
      return 0;
    }

    final File[] unmappedFiles = dir.listFiles(new SamBamGzUnmappedFileFilter());
    if (unmappedFiles == null) {
      throw new IOException("Invalid directory.");
    }

    File matedFile = null;
    File unmatedFile = null;

    for (File file : files) {
      if (file.getName().contains("unmated")) {
        if (unmatedFile != null) {
          Diagnostic.error("Multiple unmated files detected, aborting.");
          return 1;
        } else {
          unmatedFile = file;
        }
      } else {
        if (matedFile != null) {
          Diagnostic.error("Multiple mated files detected, aborting.");
          return 1;
        } else {
          matedFile = file;
        }
      }
    }

    final File unmappedFile;
    if (unmappedFiles.length > 1) {
      Diagnostic.error("Multiple unmapped files found, aborting.");
      return 1;
    } else if (unmappedFiles.length == 1) {
      unmappedFile = unmappedFiles[0];
    } else {
      unmappedFile = null;
    }

    if (matedFile == null) {
      Diagnostic.error("Unmapped file processing requires a mated file.");
      return 1;
    }
    if (unmatedFile == null) {
      Diagnostic.error("Unmapped file processing requires an unmated file.");
      return 1;
    }

    try (OutputStream rgOut = FileUtils.createOutputStream(new File(dir, UnmatedAugmenter.DEFAULT_RGSTATS_FILENAME), false)) {
      performAugmentation(matedFile, unmatedFile, unmappedFile, !mFlags.isSet(KEEP_ORIG_FLAG), true, (String) mFlags.getValue(OUTPUT_SUFFIX_FLAG), rgOut);
    }
    return 0;
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
  static void performAugmentation(File matedFile, File unmatedFile, File unmappedFile, boolean rename, boolean force, String suffix, OutputStream readGroupStatsOut) throws IOException {
    final ReadGroupStatsCalculator calc = readGroupStatsOut == null ? null : new ReadGroupStatsCalculator();
    if (matedFile != null && calc != null) {
      UnmatedAugmenter.populateReadGroupStats(matedFile, calc);
    }
    if (unmatedFile != null) {
      augment(unmatedFile, unmappedFile, rename, force, suffix, calc);
    }
    if (calc != null) {
      calc.dumpReadGroups(readGroupStatsOut);
    }
  }

  private static void augment(File unmatedFile, File unmappedFile, boolean rename, boolean force, String suffix, ReadGroupStatsCalculator calc) throws IOException {
    final File output = getOutputFile(unmatedFile, rename, force, suffix);
    final UnmatedAugmenter augmenter = new UnmatedAugmenter();
    augmenter.augmentUnmated(unmatedFile, output, FileUtils.isGzipFilename(unmatedFile), calc);

    if (unmappedFile != null) {
      final File outputUnmapped = getOutputFile(unmappedFile, rename, force, suffix);
      augmenter.augmentUnmapped(unmappedFile, outputUnmapped, FileUtils.isGzipFilename(outputUnmapped), calc);
      renameAndReindex(unmappedFile, outputUnmapped, rename);
    }

    renameAndReindex(unmatedFile, output, rename);
  }

  private static void renameAndReindex(File f, File output, boolean rename) throws IOException {
    final File finalFile;
    if (rename) {
      if (!(f.canWrite() && f.delete() && output.renameTo(f))) {
        throw new NoTalkbackSlimException("Unable to rename file: \"" + output.getPath() + "\" to: \"" + f.getPath() + "\"");
      }
      finalFile = f;
    } else {
      finalFile = output;
    }

    //recreate index if it exists
    File indexFile;
    if (SamUtils.isBAMFile(f)) {
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
    if (suffix == null) {
      throw new NullPointerException("suffix cannot be null");
    }
    final String name = f.getName().replaceAll("^(.*?)(\\.[bs]am(\\.gz)?)?$", "$1" + suffix + "$2");
    File ret = new File(f.getParentFile(), name);
    if (rename && ret.exists()) {
      ret = File.createTempFile(f.getName(), ".temp", f.getParentFile());
    } else if (ret.exists() && !force) {
      throw new NoTalkbackSlimException("File: \"" + ret.getPath() + "\" already exists. Please remove or choose a different suffix.");
    }
    return ret;
  }

  private static class UnmatedAugmenterValidator implements Validator {

    @Override
    public boolean isValid(CFlags flags) {
      final File dir = (File) flags.getAnonymousValue(0);
      if (!dir.isDirectory()) {
        flags.setParseMessage("A directory name was expected, but \"" + dir.getPath() + "\" is not a directory.");
        return false;
      }
      if (flags.isSet(OUTPUT_SUFFIX_FLAG) && !flags.isSet(KEEP_ORIG_FLAG)) {
        flags.setParseMessage("Cannot set: --" + OUTPUT_SUFFIX_FLAG + " without: --" + KEEP_ORIG_FLAG);
        return false;
      }
      return true;
    }
  }
}
